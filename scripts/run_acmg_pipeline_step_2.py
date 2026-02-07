#!/usr/bin/env python3
"""
ACMG PS1 / PM5 refinement based on Bayesian ACMG score.

This script refines PS1 and PM5 evidence by identifying amino-acid–level
recurrence among pathogenic variants (Bayesian score >= threshold)
within the same gene.

PS1 dominates PM5 when both conditions are met.
"""

import argparse
import csv
import re

import pandas as pd
from tqdm import tqdm


# ======================================================
# REGEX FOR HGVS PROTEIN NOTATION (e.g. p.Arg123Gly)
# ======================================================
AA_PATTERN = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")


# ======================================================
# CORE LOGIC
# ======================================================
def refine_ps1_pm5(
    intervar_path,
    annotation_excel,
    bayes_table_path,
    output_path,
    pathogenic_threshold
):
    # --------------------------------------------------
    # Load previous InterVar adjustments
    # --------------------------------------------------
    previous_df = pd.read_csv(intervar_path, sep="\t")

    evidence_dict = {
        f"{row.CHROM}-{row.POS}-{row.REF}-{row.ALT}": row.EVIDENCE
        for _, row in previous_df.iterrows()
    }

    # --------------------------------------------------
    # Load annotated variants (Excel)
    # --------------------------------------------------
    df = pd.read_excel(annotation_excel, index_col=0)

    df[["CHROM", "POS", "REF", "ALT"]] = (
        df[["CHROM", "POS", "REF", "ALT"]].astype(str)
    )

    # --------------------------------------------------
    # Load Bayesian score table
    # --------------------------------------------------
    df_bayes = pd.read_csv(bayes_table_path, sep="\t")

    df_bayes = (
        df_bayes[["CHROM", "POS", "REF", "ALT", "ACMG_BAYESIAN_SCORE"]]
        .drop_duplicates()
        .rename(columns={"ACMG_BAYESIAN_SCORE": "ACMG_BAYESIAN_SCORE_NEW"})
    )

    df_bayes[["CHROM", "POS", "REF", "ALT"]] = (
        df_bayes[["CHROM", "POS", "REF", "ALT"]].astype(str)
    )

    # Merge Bayesian score into main dataframe
    df = pd.merge(
        df,
        df_bayes,
        on=["CHROM", "POS", "REF", "ALT"],
        how="left"
    )

    # --------------------------------------------------
    # Missense variants only
    # --------------------------------------------------
    df_missense = df[
        df.SNPEFF_EFFECT.str.contains("missense_variant", na=False)
    ].copy()

    # --------------------------------------------------
    # PS1 / PM5 evaluation
    # --------------------------------------------------
    processed = set()

    for _, row in tqdm(
        df_missense.iterrows(),
        total=len(df_missense),
        desc="Evaluating PS1 / PM5"
    ):
        variant_key = f"{row.CHROM}-{row.POS}-{row.REF}-{row.ALT}"

        if variant_key in processed:
            continue

        if variant_key not in evidence_dict:
            continue

        processed.add(variant_key)

        if pd.isna(row.SNPEFF_AMINO_ACID_CHANGE):
            continue

        match = AA_PATTERN.match(row.SNPEFF_AMINO_ACID_CHANGE)
        if not match:
            continue

        ref_aa, position, alt_aa = match.groups()

        # Pathogenic variants in the same gene
        df_gene_pathogenic = df_missense[
            (df_missense.SNPEFF_GENE_NAME == row.SNPEFF_GENE_NAME) &
            (df_missense.ACMG_BAYESIAN_SCORE_NEW >= pathogenic_threshold)
        ]

        has_ps1 = False
        has_pm5 = False

        for _, r in df_gene_pathogenic.iterrows():

            ref_key = f"{r.CHROM}-{r.POS}-{r.REF}-{r.ALT}"
            if ref_key == variant_key:
                continue

            if pd.isna(r.SNPEFF_AMINO_ACID_CHANGE):
                continue

            r_match = AA_PATTERN.match(r.SNPEFF_AMINO_ACID_CHANGE)
            if not r_match:
                continue

            r_ref, r_pos, r_alt = r_match.groups()

            if ref_aa == r_ref and position == r_pos:
                if alt_aa == r_alt:
                    has_ps1 = True
                else:
                    has_pm5 = True

        # --------------------------------------------------
        # Final decision (PS1 overrides PM5)
        # --------------------------------------------------
        if has_ps1:
            if "PS1=1" not in evidence_dict[variant_key]:
                evidence_dict[variant_key] += ";PS1=1"
            evidence_dict[variant_key] = evidence_dict[variant_key].replace(";PM5=1", "")

        elif has_pm5:
            if "PM5=1" not in evidence_dict[variant_key]:
                evidence_dict[variant_key] += ";PM5=1"

    # --------------------------------------------------
    # Write output TSV
    # --------------------------------------------------
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["CHROM", "POS", "REF", "ALT", "EVIDENCE"])

        for key, ev in evidence_dict.items():
            chrom, pos, ref, alt = key.split("-")
            writer.writerow([chrom, pos, ref, alt, ev])

    print("✅ PS1 / PM5 refinement completed successfully.")
    print(f"📄 Output file: {output_path}")


# ======================================================
# CLI
# ======================================================
def parse_args():
    parser = argparse.ArgumentParser(
        description="Refine PS1 / PM5 using Bayesian ACMG score"
    )

    parser.add_argument(
        "--intervar",
        required=True,
        help="Input intervar_adjustments.tsv"
    )

    parser.add_argument(
        "--excel",
        required=True,
        help="Annotated Excel file (SNPEff)"
    )

    parser.add_argument(
        "--table-bayes",
        required=True,
        help="TSV table with ACMG_BAYESIAN_SCORE"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output intervar_adjustments.tsv"
    )

    parser.add_argument(
        "--pathogenic-threshold",
        type=float,
        default=0.9,
        help="Bayesian pathogenicity threshold (default: 0.9)"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    refine_ps1_pm5(
        intervar_path=args.intervar,
        annotation_excel=args.excel,
        bayes_table_path=args.table_bayes,
        output_path=args.output,
        pathogenic_threshold=args.pathogenic_threshold
    )

