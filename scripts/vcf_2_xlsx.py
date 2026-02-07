#!/usr/bin/env python3

import argparse
import pandas as pd
import sys


# ============================================================
# ANN fields according to SnpEff VCF specification
# ============================================================

ANN_FIELDS = [
    "ALLELE",
    "EFFECT",
    "IMPACT",
    "GENE",
    "GENE_ID",
    "FEATURE",
    "FEATUREID",
    "TRANSCRIPT_BIOTYPE",
    "RANK",
    "HGVS_C",
    "HGVS_P",
    "CDNA_POS",
    "CDS_POS",
    "AA_POS",
    "DISTANCE",
    "ERRORS",
]


# ============================================================
# VCF PARSING
# ============================================================

def parse_vcf(vcf_file):
    """
    Parse an annotated VCF file and expand SnpEff ANN entries
    so that each transcript annotation corresponds to one row.

    Parameters
    ----------
    vcf_file : str
        Path to the annotated VCF file.

    Returns
    -------
    pandas.DataFrame
        Expanded annotation table (one row per variant/transcript).
    """
    rows = []

    try:
        fh = open(vcf_file)
    except OSError as e:
        sys.exit(f"❌ Error opening VCF file: {e}")

    with fh:
        for line in fh:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            if len(fields) < 8:
                continue

            chrom, pos, _id, ref, alt, qual, flt, info = fields[:8]

            # ----------------------------------------------------
            # Parse INFO field
            # ----------------------------------------------------
            info_dict = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_dict[k] = v
                else:
                    info_dict[item] = True

            ann_raw = info_dict.get("ANN")
            if not ann_raw:
                continue

            ann_entries = ann_raw.split(",")

            # ----------------------------------------------------
            # Expand ANN entries (one row per transcript)
            # ----------------------------------------------------
            for ann in ann_entries:
                ann_values = ann.split("|")

                # Pad ANN values to expected length if needed
                if len(ann_values) < len(ANN_FIELDS):
                    ann_values += [""] * (len(ANN_FIELDS) - len(ann_values))

                ann_dict = dict(zip(ANN_FIELDS, ann_values))

                row = {
                    # ---- core variant fields ----
                    "CHROM": chrom,
                    "POS": int(pos),
                    "REF": ref,
                    "ALT": alt,

                    # ---- external annotations ----
                    "CLNSIG": info_dict.get("CLNSIG"),
                    "CLNREVSTAT": info_dict.get("CLNREVSTAT"),
                    "InterVarEvidence": info_dict.get("InterVarEvidence"),
                    "InterVarVeredict": info_dict.get("InterVarVeredict"),
                    "REVEL_score": info_dict.get("REVEL_score"),
                    "gnomad_exome_AF": info_dict.get("gnomad_exome_AF"),
                    "gnomad_genome_AF": info_dict.get("gnomad_genome_AF"),
                    "gnomad40_genome_faf95": info_dict.get("gnomad40_genome_faf95"),
                    "phastCons100way_vertebrate": info_dict.get("phastCons100way_vertebrate"),
                    "phastCons100way_vertebrate_rankscore": info_dict.get(
                        "phastCons100way_vertebrate_rankscore"
                    ),
                    "SAI_DS_AG": info_dict.get("SAI_DS_AG"),
                    "SAI_DS_AL": info_dict.get("SAI_DS_AL"),
                    "SAI_DS_DG": info_dict.get("SAI_DS_DG"),
                    "SAI_DS_DL": info_dict.get("SAI_DS_DL"),
                    "ACMG_BAYESIAN_SCORE": info_dict.get("ACMG_BAYESIAN_SCORE"),
                    "ACMG_BAYESIAN_VEREDICT": info_dict.get("ACMG_BAYESIAN_VEREDICT"),

                    # ---- raw ANN entry ----
                    "ANN": ann,
                }

                # ---- transcript-level fields ----
                row["ANN[*].GENE"] = ann_dict.get("GENE")
                row["ANN[*].IMPACT"] = ann_dict.get("IMPACT")
                row["ANN[*].EFFECT"] = ann_dict.get("EFFECT")
                row["ANN[*].FEATURE"] = ann_dict.get("FEATURE")
                row["ANN[*].FEATUREID"] = ann_dict.get("FEATUREID")
                row["ANN[*].RANK"] = ann_dict.get("RANK")
                row["ANN[*].HGVS_C"] = ann_dict.get("HGVS_C")
                row["ANN[*].HGVS_P"] = ann_dict.get("HGVS_P")

                rows.append(row)

    return pd.DataFrame(rows)


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Convert a SnpEff-annotated VCF into an Excel table expanded "
            "by transcript (one row per ANN entry)."
        )
    )

    parser.add_argument(
        "--vcf",
        required=True,
        help="Input VCF file annotated with SnpEff",
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output Excel (.xlsx) file",
    )

    args = parser.parse_args()

    print("📖 Reading VCF...")
    df = parse_vcf(args.vcf)

    print(f"🧬 Variant × transcript rows: {len(df)}")

    print("💾 Writing Excel file...")
    df.to_excel(args.out, index=False)

    print(f"✅ Done: {args.out}")


if __name__ == "__main__":
    main()

