#!/usr/bin/env python3

"""
Bayesian refinement of ACMG/ClinGen variant classification.

This script reads an annotated VCF containing InterVar evidence,
applies the Bayesian ACMG model (Tavtigian et al.),
and annotates each variant with posterior pathogenicity probability
and Bayesian verdict.
"""

import argparse
import pandas as pd
import subprocess


# ============================================================
# Argument parsing
# ============================================================

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Bayesian scoring of ACMG/ClinGen criteria"
    )
    parser.add_argument(
        "-v", "--vcf", required=True, help="Input annotated VCF"
    )
    parser.add_argument(
        "-o", "--out", required=True, help="Output VCF with Bayesian annotations"
    )
    parser.add_argument(
        "-e", "--evidence",
        help="Optional TSV with adjusted ACMG evidence",
        default=None
    )
    return parser.parse_args()


# ============================================================
# VCF loading utilities
# ============================================================

def load_vcf(vcf_path, parse_info=False):
    # Read VCF skipping meta-information lines
    df = pd.read_csv(
        vcf_path,
        sep="\t",
        comment="#",
        header=None,
        dtype=str
    )

    # Extract header line (#CHROM ...)
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                columns = line.lstrip("#").strip().split("\t")
                break

    df.columns = columns

    # Reload full header (## lines)
    header = subprocess.check_output(
        f"grep '^##' {vcf_path}", shell=True
    ).decode("utf-8")

    if parse_info:
        df["INFO2"] = df["INFO"].str.split(";").apply(
            lambda x: dict(
                item.split("=", 1) if "=" in item else (item, True)
                for item in x
            )
        )

    return df, header



# ============================================================
# Bayesian ACMG model
# ============================================================

PRIOR_PROBABILITY = 0.1
ODDS_PATHOGENIC_VERY_STRONG = 350


def count_evidence(row, evidence_dict):
    """Count ACMG evidence categories for a variant."""
    evidence_counter = {
        "PVS": 0, "PS": 0, "PM": 0, "PP": 0,
        "BVS": 0, "BS": 0, "BM": 0, "BP": 0,
        "BA1": 0
    }

    variant_evidence = []
    if "InterVarEvidence" in row.INFO2:
        variant_evidence = row.INFO2["InterVarEvidence"].split(",")

    variant_key = f"{row.CHROM}-{row.POS}-{row.REF}-{row.ALT}"

    # Apply external evidence adjustments
    if variant_key in evidence_dict:
        for ev in evidence_dict[variant_key]:
            label, decision = ev.split("=")
            if decision == "1" and label not in variant_evidence:
                variant_evidence.append(label)
            elif decision == "0" and label in variant_evidence:
                variant_evidence.remove(label)

    for ev in variant_evidence:
        for cat in evidence_counter:
            if ev.startswith(cat):
                evidence_counter[cat] += 1

    return evidence_counter


def calculate_posterior_probability(evidence):
    """Compute posterior probability using Tavtigian Bayesian model."""
    if evidence["BA1"] > 0:
        return 0.0

    odds = ODDS_PATHOGENIC_VERY_STRONG ** (
        evidence["PP"] / 8 +
        evidence["PM"] / 4 +
        evidence["PS"] / 2 +
        evidence["PVS"] -
        evidence["BVS"] -
        evidence["BS"] / 2 -
        evidence["BM"] / 4 -
        evidence["BP"] / 8
    )

    posterior = (odds * PRIOR_PROBABILITY) / (
        (odds - 1) * PRIOR_PROBABILITY + 1
    )

    return round(posterior, 3)


def bayesian_verdict(probability):
    """Assign ACMG classification from posterior probability."""
    if probability < 0.001:
        return "Benign"
    elif probability < 0.1:
        return "Likely_benign"
    elif probability < 0.9:
        return "Uncertain_significance"
    elif probability < 0.99:
        return "Likely_pathogenic"
    else:
        return "Pathogenic"


# ============================================================
# Main
# ============================================================

def main():
    args = parse_arguments()

    vcf_df, header = load_vcf(args.vcf, parse_info=True)
    final_vcf, _ = load_vcf(args.vcf, parse_info=False)

    vcf_df.rename(columns={"#CHROM": "CHROM"}, inplace=True)

    evidence_dict = {}
    if args.evidence:
        evidence_df = pd.read_csv(args.evidence, sep="\t")
        for _, row in evidence_df.iterrows():
            key = f"{row.CHROM}-{row.POS}-{row.REF}-{row.ALT}"
            evidence_dict[key] = row.EVIDENCE.split(";")

    print("> Calculating Bayesian ACMG scores")

    for idx, row in vcf_df.iterrows():
        if "InterVarVeredict" not in row.INFO2:
            continue

        evidence = count_evidence(row, evidence_dict)
        prob = calculate_posterior_probability(evidence)
        verdict = bayesian_verdict(prob)

        final_vcf.at[idx, "INFO"] += (
            f";ACMG_BAYESIAN_SCORE={prob}"
            f";ACMG_BAYESIAN_VEREDICT={verdict}"
        )

    with open(args.out, "w") as out:
        out.write(header)
        out.write(
            "##INFO=<ID=ACMG_BAYESIAN_SCORE,Number=1,Type=Float,"
            "Description=\"Bayesian ACMG posterior probability\">\n"
        )
        out.write(
            "##INFO=<ID=ACMG_BAYESIAN_VEREDICT,Number=1,Type=String,"
            "Description=\"Bayesian ACMG classification\">\n"
        )
        out.write("\t".join(final_vcf.columns) + "\n")

        for _, row in final_vcf.iterrows():
            out.write("\t".join(row.astype(str)) + "\n")

    print("✓ Bayesian annotation completed")


if __name__ == "__main__":
    main()
