#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Using SWOffinder outputs (<prefix>_offtargets_*.csv) and
a gene model GFF3,

- Annotate genomic regions corresponding to off-target coordinates
  (CDS / UTR / exon / intron / promoter / intergenic)
- Total edit count (#Edit)
- Number of mismatches in the seed region (10 nt immediately upstream of PAM)
- Risk category (High / Medium / Low)

The script outputs "_annotated.csv" files.

Multiple CSV files are processed in parallel.
"""

import os
import glob
import csv
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Tuple, Any

import argparse
import pandas as pd


# ===== GFF3 parsing functions =====

def parse_attributes(attr_str: str) -> Dict[str, str]:
    """Convert GFF3 attributes column into a dict"""
    d = {}
    for part in attr_str.split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d


def parse_gff3(
    gff_path: str,
) -> Tuple[Dict[str, List[dict]], Dict[str, List[dict]]]:
    """
    Parse a GFF3 file and return:
    - genes_by_chrom: {chrom: [ {start, end, strand, gene_id, attrs}, ... ]}
    - features_by_chrom: {chrom: [ {start, end, strand, type, attrs}, ... ]}
      (type: exon, CDS, five_prime_UTR, three_prime_UTR)
    """
    genes_by_chrom: Dict[str, List[dict]] = {}
    features_by_chrom: Dict[str, List[dict]] = {}

    with open(gff_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue

            attr_dict = parse_attributes(attrs)

            # gene features
            if ftype == "gene":
                gene_id = attr_dict.get("ID", attr_dict.get("Name", "NA"))
                genes_by_chrom.setdefault(seqid, []).append(
                    {
                        "start": start_i,
                        "end": end_i,
                        "strand": strand,
                        "gene_id": gene_id,
                        "attrs": attr_dict,
                    }
                )

            # exon / CDS / UTR features
            if ftype in ("exon", "CDS", "five_prime_UTR", "three_prime_UTR"):
                features_by_chrom.setdefault(seqid, []).append(
                    {
                        "start": start_i,
                        "end": end_i,
                        "strand": strand,
                        "type": ftype,
                        "attrs": attr_dict,
                    }
                )

    # Sort by start coordinate (for early break in later scans)
    for chrom in genes_by_chrom:
        genes_by_chrom[chrom].sort(key=lambda x: x["start"])
    for chrom in features_by_chrom:
        features_by_chrom[chrom].sort(key=lambda x: x["start"])

    return genes_by_chrom, features_by_chrom


# ===== Count seed mismatches from alignment =====

def count_seed_mismatches(
    aligned_target: str,
    aligned_text: str,
    seed_len: int = 10,
    pam_len: int = 3,
) -> int:
    """
    Using SWOffinder outputs:
    - aligned_target: sgRNA+PAM side
    - aligned_text: genomic side

    Count the number of mismatches in the seed_len bases
    immediately upstream of PAM.

    Roughly:
    - Skip pam_len characters from the right (assumed PAM)
    - For the seed_len characters before that, count a!=b
      (including '-' as mismatch)
    """
    if len(aligned_target) != len(aligned_text):
        # As a fallback when lengths disagree, return 0
        return 0

    i = len(aligned_target) - 1
    pam_skipped = 0
    # Skip PAM
    while i >= 0 and pam_skipped < pam_len:
        i -= 1
        pam_skipped += 1

    mismatches = 0
    seed_checked = 0
    while i >= 0 and seed_checked < seed_len:
        a = aligned_target[i]
        b = aligned_text[i]
        if a != b:
            # Treat '-' and 'N' as mismatches as well
            mismatches += 1
        seed_checked += 1
        i -= 1

    return mismatches


# ===== Genomic region annotation =====

def interval_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """Return True if [a_start, a_end] and [b_start, b_end] overlap"""
    return not (a_end < b_start or a_start > b_end)


def annotate_region(
    chrom: str,
    start: int,
    end: int,
    strand: str,
    genes_by_chrom: Dict[str, List[dict]],
    features_by_chrom: Dict[str, List[dict]],
    promoter_upstream: int = 2000,
) -> Dict[str, Any]:
    """
    Determine whether the off-target interval [start, end] is located in:
    - CDS / 5UTR / 3UTR / exon / intron / promoter / intergenic
    and return the associated gene_id(s).
    """
    region_type = "Intergenic"
    gene_ids = set()

    genes = genes_by_chrom.get(chrom, [])
    feats = features_by_chrom.get(chrom, [])

    # 1. Overlap with genes (used later for intron calls)
    overlapping_genes = []
    for g in genes:
        if g["end"] < start:
            continue
        if g["start"] > end:
            break
        if interval_overlap(start, end, g["start"], g["end"]):
            overlapping_genes.append(g)
            gene_ids.add(g["gene_id"])

    # 2. Overlap with exon / CDS / UTR
    overlapping_types = set()
    for ft in feats:
        if ft["end"] < start:
            continue
        if ft["start"] > end:
            break
        if interval_overlap(start, end, ft["start"], ft["end"]):
            overlapping_types.add(ft["type"])

    # Decide region_type by priority
    # (Feel free to adjust)
    if "CDS" in overlapping_types:
        region_type = "CDS"
    elif "five_prime_UTR" in overlapping_types:
        region_type = "5UTR"
    elif "three_prime_UTR" in overlapping_types:
        region_type = "3UTR"
    elif "exon" in overlapping_types:
        region_type = "Exon"
    elif overlapping_genes:
        region_type = "Intron"

    # 3. If not inside any gene, check promoter / intergenic
    if not overlapping_genes:
        for g in genes:
            # Define promoter as flanking the TSS
            if g["strand"] == "+":
                prom_start = g["start"] - promoter_upstream
                prom_end = g["start"] - 1
            else:
                prom_start = g["end"] + 1
                prom_end = g["end"] + promoter_upstream
            if interval_overlap(start, end, prom_start, prom_end):
                region_type = "Promoter"
                gene_ids.add(g["gene_id"])
                overlapping_genes.append(g)

    gene_id_str = ";".join(sorted(gene_ids)) if gene_ids else ""

    return {
        "RegionType": region_type,
        "GeneIDs": gene_id_str,
    }


# ===== Function to process one CSV =====

def process_one_csv(
    csv_path: str,
    gff_path: str,
    promoter_upstream: int = 2000,
    seed_len: int = 10,
) -> str:
    """
    Read one SWOffinder CSV,
    add:
      - genomic region annotations for each off-target
      - seed mismatches
      - risk category
    and write out an annotated CSV.
    """
    print(f"[INFO] Processing {os.path.basename(csv_path)}")

    # Read GFF3 (for the small number of parallel processes,
    # re-reading in each worker is acceptable)
    genes_by_chrom, features_by_chrom = parse_gff3(gff_path)

    df = pd.read_csv(csv_path)

    # Compute genomic interval from AlignedText
    # (length on genome side, excluding '-')
    genome_lengths = df["AlignedText"].apply(
        lambda s: sum(1 for c in str(s) if c != "-")
    )
    df["GenomicLength"] = genome_lengths

    # Compute start / end regardless of strand
    starts = []
    ends = []
    for strand, end_pos, glen in zip(
        df["Strand"], df["EndPosition"], df["GenomicLength"]
    ):
        end_i = int(end_pos)
        start_i = end_i - int(glen) + 1
        starts.append(start_i)
        ends.append(end_i)

    df["GenomicStart"] = starts
    df["GenomicEnd"] = ends

    # Count seed mismatches
    seed_mismatches = []
    for at, tx in zip(df["AlignedTarget"], df["AlignedText"]):
        seed_mismatches.append(
            count_seed_mismatches(str(at), str(tx), seed_len=seed_len, pam_len=3)
        )
    df["SeedMismatches"] = seed_mismatches

    # Annotate genomic regions
    region_types = []
    gene_ids = []
    for chrom, start, end, strand in zip(
        df["Chromosome"], df["GenomicStart"], df["GenomicEnd"], df["Strand"]
    ):
        ann = annotate_region(
            chrom=str(chrom),
            start=int(start),
            end=int(end),
            strand=str(strand),
            genes_by_chrom=genes_by_chrom,
            features_by_chrom=features_by_chrom,
            promoter_upstream=promoter_upstream,
        )
        region_types.append(ann["RegionType"])
        gene_ids.append(ann["GeneIDs"])

    df["RegionType"] = region_types
    df["GeneIDs"] = gene_ids

    # Assign risk categories
    risk_levels = []
    for total_edit, seed_mm in zip(df["#Edit"], df["SeedMismatches"]):
        te = int(total_edit)
        sm = int(seed_mm)
        if te <= 2 and sm <= 1:
            risk = "High"    # high off-target risk
        elif te <= 3 and sm <= 2:
            risk = "Medium"  # moderate
        else:
            risk = "Low"     # relatively low
        risk_levels.append(risk)

    df["RiskCategory"] = risk_levels

    # Output
    out_path = os.path.splitext(csv_path)[0] + "_annotated.csv"
    df.to_csv(out_path, index=False)
    print(f"[INFO] Wrote {os.path.basename(out_path)}")
    return out_path


# ===== Argument parser =====

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate SWOffinder off-target CSV files using a GFF3 gene model."
    )
    parser.add_argument(
        "--base-dir",
        default=".",
        help="Directory containing SWOffinder CSV files (default: current directory).",
    )
    parser.add_argument(
        "--gff",
        required=True,
        help="Path to GFF3 gene model.",
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help="Prefix used by SWOffinder in output filenames "
             "(e.g. 'Dnig' for files like 'Dnig_offtargets_*.csv').",
    )
    parser.add_argument(
        "--promoter-upstream",
        type=int,
        default=2000,
        help="Upstream distance (bp) to define promoter region (default: 2000).",
    )
    parser.add_argument(
        "--seed-len",
        type=int,
        default=10,
        help="Seed length (nt) immediately upstream of PAM (default: 10).",
    )
    return parser.parse_args()


# ===== Main =====

def main():
    args = parse_args()

    base_dir = args.base_dir
    gff_path = args.gff
    promoter_upstream = args.promoter_upstream
    seed_len = args.seed_len

    pattern = os.path.join(base_dir, f"{args.prefix}_offtargets_*.csv")

    csv_files = sorted(glob.glob(pattern))
    if not csv_files:
        print("[ERROR] No CSV files found for pattern:", pattern)
        return

    print("[INFO] Base directory:", base_dir)
    print("[INFO] GFF3 file:", gff_path)
    print("[INFO] File pattern:", pattern)
    print("[INFO] Found CSV files:")
    for p in csv_files:
        print("  -", os.path.basename(p))

    # Parallel processing
    with ProcessPoolExecutor() as ex:
        futures = []
        for csv_path in csv_files:
            futures.append(
                ex.submit(
                    process_one_csv,
                    csv_path,
                    gff_path,
                    promoter_upstream,  # upstream distance for promoter definition
                    seed_len,           # seed length (nt immediately upstream of PAM)
                )
            )
        for fut in futures:
            try:
                fut.result()
            except Exception as e:
                print("[ERROR] worker failed:", e)

    print("[INFO] Done.")


if __name__ == "__main__":
    main()
