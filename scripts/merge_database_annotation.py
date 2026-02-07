#!/usr/bin/env python3

"""
Merge expanded VCF annotations with a curated variant database
using canonical (MANE) transcripts.

This script:
- Filters transcript-expanded annotations to canonical transcripts
- Restricts annotations to genes present in the curated database
- Merges annotations back into the curated variant table
- Optionally exports the filtered annotation table for QC/debugging
"""

import argparse
import json
import pandas as pd


# ============================================================
# Argument parsing
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge curated variant database with transcript-aware annotations"
    )

    parser.add_argument(
        "-b", "--base",
        required=True,
        help="Curated variant Excel file (e.g. rs_base.xlsx)"
    )

    parser.add_argument(
        "-a", "--annotation",
        required=True,
        help="Transcript-expanded annotation Excel file"
    )

    parser.add_argument(
        "-t", "--transcripts",
        required=True,
        help="JSON file with canonical (MANE) transcripts per gene"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output merged Excel file"
    )

    parser.add_argument(
        "--annot-filtered",
        default=None,
        help="Optional path to save filtered annotation table (QC file)"
    )

    return parser.parse_args()


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    # --------------------------------------------------------
    # Load curated database
    # --------------------------------------------------------

    df_base = pd.read_excel(args.base)

    print(f"Curated database loaded: {df_base.shape[0]} rows")

    # --------------------------------------------------------
    # Load canonical transcripts
    # --------------------------------------------------------

    with open(args.transcripts, "r") as fh:
        transcripts_dict = json.load(fh)

    # Remove version suffix from transcript IDs
    transcripts_clean = {
        gene: transcript.split(".")[0]
        for gene, transcript in transcripts_dict.items()
    }

    # --------------------------------------------------------
    # Load annotation table (only required columns)
    # --------------------------------------------------------

    usecols = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "CLNSIG",
        "CLNREVSTAT",
        "InterVarEvidence",
        "InterVarVeredict",
        "REVEL_score",
        "gnomad_exome_AF",
        "gnomad_genome_AF",
        "gnomad40_genome_faf95",
        "SAI_DS_AG",
        "SAI_DS_AL",
        "SAI_DS_DG",
        "SAI_DS_DL",
        "ACMG_BAYESIAN_SCORE",
        "ACMG_BAYESIAN_VEREDICT",
        "ANN",
        "ANN[*].GENE",
        "ANN[*].IMPACT",
        "ANN[*].EFFECT",
        "ANN[*].FEATURE",
        "ANN[*].FEATUREID",
        "ANN[*].RANK",
        "ANN[*].HGVS_C",
        "ANN[*].HGVS_P",
    ]

    df_annot = pd.read_excel(args.annotation, usecols=usecols)

    print(f"Annotation table loaded: {df_annot.shape[0]} rows")

    # --------------------------------------------------------
    # Normalize transcript IDs
    # --------------------------------------------------------

    df_annot["TRANSCRIPT_CLEAN"] = (
        df_annot["ANN[*].FEATUREID"]
        .astype(str)
        .str.split(".")
        .str[0]
    )

    # --------------------------------------------------------
    # Filter to canonical transcripts
    # --------------------------------------------------------
    # IMPORTANT:
    # Variants without a canonical transcript here WILL BE LOST
    # --------------------------------------------------------

    canonical_mask = df_annot.apply(
        lambda row: (
            row["ANN[*].GENE"] in transcripts_clean
            and row["TRANSCRIPT_CLEAN"] == transcripts_clean[row["ANN[*].GENE"]]
        ),
        axis=1
    )

    df_annot_filt = df_annot[canonical_mask].copy()

    # --------------------------------------------------------
    # Restrict to genes present in curated database
    # --------------------------------------------------------

    curated_genes = set(df_base["gen"].dropna())

    df_annot_filt = df_annot_filt[
        df_annot_filt["ANN[*].GENE"].isin(curated_genes)
    ].copy()

    # --------------------------------------------------------
    # Remove duplicate variants
    # --------------------------------------------------------
    # NOTE:
    # This collapses multiple transcripts per variant
    # --------------------------------------------------------

    df_annot_filt.drop_duplicates(
        subset=["CHROM", "POS", "REF", "ALT"],
        inplace=True
    )

    print(f"Filtered annotation table: {df_annot_filt.shape[0]} rows")

    # --------------------------------------------------------
    # Save filtered annotation (optional QC output)
    # --------------------------------------------------------

    if args.annot_filtered:
        df_annot_filt.to_excel(args.annot_filtered, index=False)
        print(f"Filtered annotation saved to: {args.annot_filtered}")

    # --------------------------------------------------------
    # Harmonize merge key types
    # --------------------------------------------------------

    for col in ["CHROM", "REF", "ALT"]:
        df_annot_filt[col] = df_annot_filt[col].astype(str)

    df_annot_filt["POS"] = df_annot_filt["POS"].astype(int)

    for col in ["Chr", "Ref", "Alt"]:
        df_base[col] = df_base[col].astype(str)

    df_base["Pos"] = df_base["Pos"].astype(int)

    # --------------------------------------------------------
    # Merge curated database with annotations
    # --------------------------------------------------------

    df_merged = df_base.merge(
        df_annot_filt,
        how="left",
        left_on=["Chr", "Pos", "Ref", "Alt"],
        right_on=["CHROM", "POS", "REF", "ALT"],
        suffixes=("", "_ann")
    )

    # --------------------------------------------------------
    # Create SnpEff-style columns
    # --------------------------------------------------------

    df_merged["SNPEFF_AMINO_ACID_CHANGE"] = df_merged["ANN[*].HGVS_P"]
    df_merged["SNPEFF_CODON_CHANGE"] = df_merged["ANN[*].HGVS_C"]
    df_merged["SNPEFF_EFFECT"] = df_merged["ANN[*].EFFECT"]
    df_merged["SNPEFF_GENE_NAME"] = df_merged["ANN[*].GENE"]
    df_merged["SNPEFF_IMPACT"] = df_merged["ANN[*].IMPACT"]
    df_merged["SNPEFF_TRANSCRIPT_ID"] = df_merged["ANN[*].FEATUREID"]

    df_merged["TRANSCRIPT_CLEAN"] = (
        df_merged["SNPEFF_TRANSCRIPT_ID"]
        .astype(str)
        .str.split(".")
        .str[0]
    )

    # --------------------------------------------------------
    # Drop redundant columns
    # --------------------------------------------------------

    cols_to_drop = [
        "Chr", "Pos", "Ref", "Alt",
        "ANN[*].GENE",
        "ANN[*].IMPACT",
        "ANN[*].EFFECT",
        "ANN[*].FEATUREID",
        "ANN[*].HGVS_C",
        "ANN[*].HGVS_P",
    ]

    df_merged.drop(
        columns=[c for c in cols_to_drop if c in df_merged.columns],
        inplace=True
    )

    # --------------------------------------------------------
    # Save final merged database
    # --------------------------------------------------------

    df_merged.to_excel(args.output, index=False)

    annotated = df_merged["SNPEFF_GENE_NAME"].notna().sum()
    total = df_merged.shape[0]

    print("\n✅ Merge completed successfully")
    print(f"Output file: {args.output}")
    print(f"Annotated variants: {annotated}/{total} ({annotated/total:.1%})")


# ============================================================
if __name__ == "__main__":
    main()

