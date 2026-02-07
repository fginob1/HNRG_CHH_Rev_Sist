#!/usr/bin/env python3
"""
Merge refined Bayesian ACMG classifications into the annotated variant database.

This script integrates the final Bayesian ACMG score and verdict (from a TSV
generated with GATK VariantsToTable) into the main annotated Excel database.
"""

import pandas as pd
import argparse


def load_inputs(excel_path, bayes_table_path):
    """Load Excel variant database and Bayesian TSV table."""
    df_excel = pd.read_excel(excel_path)
    df_bayes = pd.read_csv(bayes_table_path, sep="\t")
    return df_excel, df_bayes


def normalize_variant_columns(df):
    """Ensure standard variant coordinate column names."""
    return df.rename(columns={
        "Chr": "CHROM",
        "Pos": "POS",
        "Ref": "REF",
        "Alt": "ALT"
    })


def prepare_bayes_table(df_bayes):
    """
    Keep one row per variant and rename Bayesian columns
    to indicate refined classification.
    """
    df_bayes = df_bayes.rename(columns={
        "ACMG_BAYESIAN_SCORE": "ACMG_BAYESIAN_SCORE_REFINED",
        "ACMG_BAYESIAN_VEREDICT": "ACMG_BAYESIAN_VEREDICT_REFINED"
    })

    return df_bayes.drop_duplicates(
        subset=["CHROM", "POS", "REF", "ALT"]
    )


def merge_databases(df_excel, df_bayes):
    """Left-merge refined Bayesian annotations into Excel database."""
    return df_excel.merge(
        df_bayes[
            [
                "CHROM", "POS", "REF", "ALT",
                "ACMG_BAYESIAN_SCORE_REFINED",
                "ACMG_BAYESIAN_VEREDICT_REFINED"
            ]
        ],
        on=["CHROM", "POS", "REF", "ALT"],
        how="left"
    )


def main(args):

    # ==========================
    # LOAD DATA
    # ==========================
    df_excel, df_bayes = load_inputs(args.excel, args.table_bayes)

    print(f"📊 Original Excel shape: {df_excel.shape}")
    print(f"📄 Bayesian table shape: {df_bayes.shape}")

    # ==========================
    # NORMALIZE COLUMNS
    # ==========================
    df_excel = normalize_variant_columns(df_excel)
    df_bayes = prepare_bayes_table(df_bayes)

    # ==========================
    # MERGE
    # ==========================
    df_merged = merge_databases(df_excel, df_bayes)

    print(f"🔗 Merged dataset shape: {df_merged.shape}")

    # ==========================
    # WRITE OUTPUTS
    # ==========================
    if args.output_tsv:
        df_merged.to_csv(args.output_tsv, sep="\t", index=False)
        print(f"📄 TSV written: {args.output_tsv}")

    if args.output_xlsx:
        df_merged.to_excel(args.output_xlsx, index=False)
        print(f"📊 Excel written: {args.output_xlsx}")

    print("✅ Refinement merge completed successfully")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Merge refined Bayesian ACMG classifications into variant database"
    )

    parser.add_argument(
        "--excel",
        required=True,
        help="Annotated variant Excel database (bdd_annotated.xlsx)"
    )

    parser.add_argument(
        "--table-bayes",
        required=True,
        help="Bayesian TSV table (from GATK VariantsToTable)"
    )

    parser.add_argument(
        "--output-tsv",
        default=None,
        help="Optional output TSV file"
    )

    parser.add_argument(
        "--output-xlsx",
        default=None,
        help="Optional output Excel file"
    )

    args = parser.parse_args()
    main(args)

