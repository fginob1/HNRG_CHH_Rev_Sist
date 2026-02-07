#!/usr/bin/env python3
# ======================================================
# ACMG PM3 assignment script
# ======================================================

import pandas as pd
import csv
import argparse
from tqdm import tqdm

# ======================================================
# MAIN
# ======================================================

def main(args):

    # ==================================================
    # LOAD PREVIOUS INTERVAR ADJUSTMENTS
    # ==================================================
    previous_df = pd.read_csv(args.intervar, sep="\t")
    evidence_dict = {
        f"{r.CHROM}-{r.POS}-{r.REF}-{r.ALT}": r.EVIDENCE
        for _, r in previous_df.iterrows()
    }

    # ==================================================
    # LOAD ANNOTATED EXCEL
    # ==================================================
    df = pd.read_excel(args.excel, index_col=0)

    df["pmid"] = df["pmid"].astype(str).str.strip()
    df["id_paper"] = df["id_paper"].astype(str).str.strip()
    df[["CHROM", "POS", "REF", "ALT"]] = df[["CHROM", "POS", "REF", "ALT"]].astype(str)

    # ==================================================
    # LOAD BAYESIAN TABLE (SOURCE OF TRUTH)
    # ==================================================
    df_bayes = pd.read_csv(args.table_bayes, sep="\t")

    df_bayes = (
        df_bayes[
            [
                "CHROM", "POS", "REF", "ALT",
                "ACMG_BAYESIAN_SCORE",
                "ACMG_BAYESIAN_VEREDICT"
            ]
        ]
        .drop_duplicates()
        .rename(
            columns={
                "ACMG_BAYESIAN_SCORE": "ACMG_BAYESIAN_SCORE_NEW",
                "ACMG_BAYESIAN_VEREDICT": "ACMG_BAYESIAN_VERDICT_NEW"
            }
        )
    )

    df_bayes[["CHROM", "POS", "REF", "ALT"]] = df_bayes[
        ["CHROM", "POS", "REF", "ALT"]
    ].astype(str)

    df = pd.merge(
        df,
        df_bayes,
        on=["CHROM", "POS", "REF", "ALT"],
        how="left"
    )

    # ==================================================
    # LOAD GENE INHERITANCE TABLE
    # ==================================================
    gene_inh_df = pd.read_csv(args.gene_inh, sep="\t", index_col="Gene")

    def get_inheritance_models(gene):
        if gene not in gene_inh_df.index:
            return []
        v = gene_inh_df.loc[gene, "Inheritance_model"]
        if pd.isna(v):
            return []
        return v.split(",")

    # ==================================================
    # PATIENT / VARIANT OBJECTS
    # ==================================================
    class Variant:
        def __init__(self, row):
            self.gene = row.SNPEFF_GENE_NAME
            self.coord = f"{row.CHROM}-{row.POS}-{row.REF}-{row.ALT}"
            self.zygosity = row.cigocity.strip()
            self.faf = row.gnomad40_genome_faf95
            self.score = row.ACMG_BAYESIAN_SCORE_NEW
            self.verdict = row.ACMG_BAYESIAN_VERDICT_NEW

    class Patient:
        def __init__(self):
            self.variants = []

    # ==================================================
    # BUILD PATIENT STRUCTURE
    # ==================================================
    patients = []

    for pmid in tqdm(df.pmid.unique(), desc="Building patients"):
        df_pmid = df[df.pmid == pmid]
        for pid in df_pmid.id_paper.unique():
            patient = Patient()
            rows = df_pmid[df_pmid.id_paper == pid]
            for _, r in rows.iterrows():
                patient.variants.append(Variant(r))
            patients.append(patient)

    print(f"🧬 Patients loaded: {len(patients)}")

    # ==================================================
    # INITIALIZE PM3 COUNTERS
    # ==================================================
    variants_dict = {}

    for patient in patients:
        for v in patient.variants:
            if (
                (pd.isna(v.faf) or v.faf <= args.mcpaf)
                and "AR" in get_inheritance_models(v.gene)
            ):
                key = f"{v.gene}: {v.coord}"
                if key not in variants_dict:
                    variants_dict[key] = {
                        "Het_PLP": 0.0,
                        "Het_VUS": 0.0,
                        "Homo": 0.0
                    }

    # ==================================================
    # VARIANT FREQUENCY IN COHORT
    # ==================================================
    variant_counts = {}
    for patient in patients:
        for v in patient.variants:
            key = f"{v.gene}: {v.coord}"
            variant_counts[key] = variant_counts.get(key, 0) + 1

    # ==================================================
    # PM3 SCORING
    # ==================================================
    for key in tqdm(variants_dict, desc="Calculating PM3"):
        gene, coord = key.split(": ")

        for patient in patients:
            gene_variants = [v for v in patient.variants if v.gene == gene]
            if not gene_variants:
                continue

            # Rank variants by Bayesian score, then cohort frequency
            ranked = sorted(
                gene_variants,
                key=lambda v: (
                    -float(v.score) if not pd.isna(v.score) else 0,
                    -variant_counts.get(f"{v.gene}: {v.coord}", 0)
                )
            )[:2]

            if coord not in [v.coord for v in ranked]:
                continue

            alleles = 0

            for v in ranked:
                if alleles >= 2:
                    break

                if v.coord == coord:
                    if v.zygosity == "Homo":
                        alleles += 2
                        if variants_dict[key]["Homo"] < 1:
                            variants_dict[key]["Homo"] += 0.5
                    elif v.zygosity == "Hetero":
                        alleles += 1
                else:
                    if v.zygosity == "Hetero":
                        alleles += 1
                        if v.verdict in ["Pathogenic", "Likely_pathogenic"]:
                            variants_dict[key]["Het_PLP"] += 1
                        elif v.verdict == "Uncertain_significance":
                            if variants_dict[key]["Het_VUS"] < 0.5:
                                variants_dict[key]["Het_VUS"] += 0.25

    # ==================================================
    # ASSIGN PM3 GRADES
    # ==================================================
    for key, scores in variants_dict.items():
        coord = key.split(": ")[1]
        points = scores["Het_PLP"] + scores["Het_VUS"] + scores["Homo"]

        if coord not in evidence_dict:
            continue

        if points >= 4:
            evidence_dict[coord] += ";PM3=1;grade_PM3=0"
        elif points >= 2:
            evidence_dict[coord] += ";PM3=1;grade_PM3=1"
        elif points >= 1:
            evidence_dict[coord] += ";PM3=1;grade_PM3=2"
        elif points >= 0.5:
            evidence_dict[coord] += ";PM3=1;grade_PM3=3"
        else:
            evidence_dict[coord] += ";PM3=0"

    # ==================================================
    # WRITE OUTPUT
    # ==================================================
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["CHROM", "POS", "REF", "ALT", "EVIDENCE"])
        for k, v in evidence_dict.items():
            chrom, pos, ref, alt = k.split("-")
            writer.writerow([chrom, pos, ref, alt, v])

    print("✅ PM3 assignment completed successfully")
    print(f"📄 Output file: {args.output}")


# ======================================================
# CLI
# ======================================================
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="ACMG PM3 assignment (compound heterozygosity / homozygosity)"
    )

    parser.add_argument("--excel", required=True, help="Annotated Excel file")
    parser.add_argument("--table-bayes", required=True, help="Bayesian TSV table")
    parser.add_argument("--gene-inh", required=True, help="Gene inheritance TSV")
    parser.add_argument("--intervar", required=True, help="Previous intervar_adjustments.tsv")
    parser.add_argument("--output", required=True, help="Output intervar_adjustments.tsv")
    parser.add_argument("--mcpaf", type=float, required=True, help="Maximum credible population AF")

    args = parser.parse_args()
    main(args)

