#!/usr/bin/env python3

import pandas as pd
import numpy as np
import math
import json
import csv
import argparse
from tqdm import tqdm
from collections import Counter

# ============================================================
# CLASES
# ============================================================

class patient_object(object):
    def __init__(self, row):
        self.Patient_ID = row.row_id
        self.PMID = row.pmid
        self.ID = row.id_paper
        self.Year = row.year_of_public
        self.Consanguinity = row.consanguinity
        self.Method = row.method
        self.Sex = row.SEX.strip().upper()
        self.Variants = []

    def add_variant(self, variant):
        self.Variants.append(variant)
        self.Variants.sort(key=lambda x: x.ACMG_Bayesian_score, reverse=True)

class variant_object(object):
    def __init__(self, row):
        self.CHROM = row.CHROM
        self.POS = row.POS
        self.REF = row.REF
        self.ALT = row.ALT
        self.Zygosity = row.cigocity.strip()
        self.AA_change = "." if pd.isna(row.SNPEFF_AMINO_ACID_CHANGE) else row.SNPEFF_AMINO_ACID_CHANGE
        self.REVEL = row.REVEL_score
        self.SpliceAI_acc_gain = float(max(str(row.SAI_DS_AG).split(',')))
        self.SpliceAI_acc_loss = float(max(str(row.SAI_DS_AL).split(',')))
        self.SpliceAI_don_gain = float(max(str(row.SAI_DS_DG).split(',')))
        self.SpliceAI_don_loss = float(max(str(row.SAI_DS_DL).split(',')))
        self.FAF = row.gnomad40_genome_faf95
        self.InterVar_labels = [] if pd.isna(row.InterVarEvidence) else row.InterVarEvidence.split(',')
        self.Gene = row.SNPEFF_GENE_NAME
        self.Effect = row.SNPEFF_EFFECT.split('&')
        self.ACMG_Bayesian_score = 0.5 if pd.isna(row.ACMG_BAYESIAN_SCORE) else row.ACMG_BAYESIAN_SCORE

    @property
    def key(self):
        return f"{self.CHROM}-{self.POS}-{self.REF}-{self.ALT}"

def get_inheritance_models(gene, gene_inh_df):
    if gene not in gene_inh_df.index:
        return []
    v = gene_inh_df.loc[gene, "Inheritance_model"]
    if pd.isna(v):
        return []
    return v.split(",")

def chrom_sort_key(chrom):
    chrom = str(chrom)
    if chrom.isdigit():
        return int(chrom)
    if chrom.upper() == "X":
        return 23
    if chrom.upper() == "Y":
        return 24
    return 25

def write_intervar(evidence_dict, output_path):
    sorted_keys = sorted(
        evidence_dict.keys(),
        key=lambda k: (chrom_sort_key(k.split("-")[0]), int(k.split("-")[1]))
    )

    with open(output_path, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["CHROM", "POS", "REF", "ALT", "EVIDENCE"])
        for k in sorted_keys:
            chrom, pos, ref, alt = k.split("-")
            writer.writerow([chrom, pos, ref, alt, evidence_dict[k]])

def compute_maximum_credible_AF(
    patients,
    variants,
    gene_inh_df,
    inheritance_model,
    prevalence,
    penetrance
):

    patients_number = len(patients)

    # ==========================
    # AUTOSÓMICO RECESIVO (AR)
    # ==========================
    if inheritance_model == "AR":

        genes = sorted(set(v.Gene for v in variants))

        patients_genes = {g: 0 for g in genes}
        for p in patients:
            for g in {v.Gene for v in p.Variants}:
                patients_genes[g] += 1

        patients_genes_sorted = dict(
            sorted(patients_genes.items(), key=lambda x: x[1], reverse=True)
        )

        max_gene = None
        for gene in patients_genes_sorted:
            inh = gene_inh_df.loc[gene, "Inheritance_model"] if gene in gene_inh_df.index else ""
            if isinstance(inh, str) and "AR" in inh:
                max_gene = gene
                break

        if max_gene is None:
            raise RuntimeError("No se encontró ningún gen AR")

        max_genetic_contribution = patients_genes_sorted[max_gene] / patients_number
        genetic_se = math.sqrt(
            max_genetic_contribution * (1 - max_genetic_contribution) / patients_number
        )

        allele_counter = Counter()
        for p in patients:
            for v in p.Variants:
                if v.Gene == max_gene:
                    allele_counter[v.key] += 2 if v.Zygosity == "Homo" else 1

        max_variant = allele_counter.most_common(1)[0][0]
        total_alleles = sum(allele_counter.values())

        max_allelic_contribution = allele_counter[max_variant] / total_alleles
        allelic_se = math.sqrt(
            max_allelic_contribution * (1 - max_allelic_contribution) / patients_number
        )

        maximum_credible_population_AF = (
            math.sqrt(prevalence)
            * (max_allelic_contribution + 3 * allelic_se)
            * math.sqrt(max_genetic_contribution + 3 * genetic_se)
            / math.sqrt(penetrance)
        )

        print(f"Maximum credible population allele frequency computed: {maximum_credible_population_AF}")

        return maximum_credible_population_AF, max_gene, max_variant

    # ==========================
    # AUTOSÓMICO DOMINANTE (AD)
    # ==========================
    elif inheritance_model == "AD":

        variant_counter = Counter()
        for p in patients:
            for v in p.Variants:
                variant_counter[v.key] += 1

        max_variant = variant_counter.most_common(1)[0][0]
        max_allelic_contribution = variant_counter[max_variant] / patients_number

        std_error = math.sqrt(
            max_allelic_contribution * (1 - max_allelic_contribution) / patients_number
        )

        maximum_credible_population_AF = (
            prevalence * (max_allelic_contribution + 3 * std_error) / penetrance
        )

        print(f"Maximum credible population allele frequency computed: {maximum_credible_population_AF}")

        return maximum_credible_population_AF, None, max_variant

    else:
        raise ValueError("Modelo de herencia inválido")

# ============================================================
# PIPELINE PRINCIPAL
# ============================================================

def run_pipeline(
    path_excel_variantes,
    path_gene_inh_tsv,
    path_gnomad_constraint,
    path_transcripts_json,
    output_intervar_adjustments,
    prevalence,
    penetrance,
    inheritance_model
):

    # ---------------- Script 1: carga de pacientes ----------------
    df = pd.read_excel(path_excel_variantes, index_col=0)
    df['pmid'] = df['pmid'].astype(str).str.strip()
    df['id_paper'] = df['id_paper'].astype(str).str.strip()

    patients = []
    for pmid in tqdm(df.pmid.unique(), desc='Cargando pacientes'):
        for pid in df[df.pmid == pmid].id_paper.unique():
            rows = df[(df.pmid == pmid) & (df.id_paper == pid)]
            patient = patient_object(rows.iloc[0])
            for _, r in rows.iterrows():
                patient.add_variant(variant_object(r))
            patients.append(patient)

    # ---------------- variantes únicas ----------------
    variant_dict = {}
    for p in patients:
        for v in p.Variants:
            if v.key not in variant_dict:
                variant_dict[v.key] = v

    variants = list(variant_dict.values())

    # ---------------- Script 2: genes ----------------
    genes_list = sorted(set(v.Gene for v in variants))

    # ---------------- Script 3: herencia ----------------
    gene_inh_df = pd.read_csv(path_gene_inh_tsv, sep='\t', index_col='Gene')

    # ---------------- Script 4: maximum credible population AF (AR, conservador) ----------------
    maximum_credible_population_AF, max_gene, max_variant = compute_maximum_credible_AF(
    patients=patients,
    variants=variants,
    gene_inh_df=gene_inh_df,
    inheritance_model=inheritance_model,
    prevalence=prevalence,
    penetrance=penetrance)

    # ---------------- Script 5: inicialización PM2 ----------------
    evidence = {k: 'grade_PM2=3' for k in variant_dict.keys()}

    # ---------------- Script 6: BA1 / BS1 / PM2 ----------------
    for v in variants:
        if v.FAF >= maximum_credible_population_AF * 10:
            evidence[v.key] += ';BA1=1;BS1=0;PM2=0'
        elif v.FAF > maximum_credible_population_AF:
            evidence[v.key] += ';BA1=0;BS1=1;PM2=0'
        else:
            evidence[v.key] += ';BA1=0;BS1=0;PM2=1'

    # ---------------- Script 7: REVEL (lógica completa original) ----------------

    for v in variants:

        if "missense_variant" not in v.Effect:
            continue

        ev = evidence[v.key]
        labels = v.InterVar_labels

        # evitar duplicación
        if "PP3=" in ev or "BP4=" in ev:
            continue

        r = v.REVEL

        # --------- PP3 ---------
        if r >= 0.932:
            if "PM1" not in labels:
                evidence[v.key] += ";PP3=1;grade_PP3=1;BP4=0"
            else:
                evidence[v.key] += ";PP3=1;grade_PP3=2;BP4=0"

        elif 0.932 > r >= 0.773:
            evidence[v.key] += ";PP3=1;grade_PP3=2;BP4=0"

        elif 0.773 > r >= 0.644:
            evidence[v.key] += ";PP3=1;grade_PP3=2;BP4=0"

        elif 0.644 > r > 0.290:
            evidence[v.key] += ";PP3=0;BP4=0"

        # --------- BP4 ---------
        elif 0.290 >= r > 0.183:
            evidence[v.key] += ";BP4=1;grade_BP4=3;PP3=0"

        elif 0.183 >= r > 0.016:
            evidence[v.key] += ";BP4=1;grade_BP4=2;PP3=0"

        elif 0.016 >= r > 0.003:
            evidence[v.key] += ";BP4=1;grade_BP4=1;PP3=0"

        elif r <= 0.003:
            evidence[v.key] += ";BP4=1;grade_BP4=0;PP3=0"


    # ---------------- Script 8: PS4 ----------------
    ps4_list = []
    for p in patients:
        for v in p.Variants:
            if v.FAF <= maximum_credible_population_AF or math.isnan(v.FAF):
                inh = get_inheritance_models(v.Gene, gene_inh_df)
                if any(x in inh for x in ['AD', 'XR', 'XL']):
                    if v.Zygosity in ['Hetero', 'Hemi']:
                        ps4_list.append(f'{v.Gene}:{v.key}')

    counts = Counter(ps4_list)
    for k, c in counts.items():
        coord = k.split(':')[1]
        if c >= 16:
            evidence[coord] += ';PS4=1;grade_PS4=0'
        elif c >= 8:
            evidence[coord] += ';PS4=1;grade_PS4=1'
        elif c >= 4:
            evidence[coord] += ';PS4=1;grade_PS4=2'
        elif c >= 2:
            evidence[coord] += ';PS4=1;grade_PS4=3'
        else:
            evidence[coord] += ';PS4=0'

    # ---------------- Script 9: PP2 ----------------
    gnomad_df = pd.read_csv(path_gnomad_constraint, sep="\t")
    gnomad_df = gnomad_df[gnomad_df.transcript.str.startswith("ENST")]

    with open(path_transcripts_json) as f:
        transcripts = json.load(f)

    # cache de z-scores
    zscore_dict = dict(
        zip(gnomad_df["transcript"], gnomad_df["mis.z_score"])
    )

    for v in variants:

        ev = evidence[v.key]
        labels = v.InterVar_labels

        # solo missense
        if "missense_variant" not in v.Effect:
            continue

        # cromosoma X
        if v.CHROM == "X":
            if "PP2=" not in ev:
                evidence[v.key] += ";PP2=0"
            continue

        # transcript
        enst = transcripts.get(v.Gene, "").split(".")[0]
        gen_z = zscore_dict.get(enst, None)

        # gen no constrained
        if gen_z is None or gen_z < 3.09:
            if "PP2=" not in ev:
                evidence[v.key] += ";PP2=0"
            continue

        # -------- gen constrained --------

        # PP3 = 0
        if "PP3=0" in ev:
            if "PP2=" not in ev:
                evidence[v.key] += ";PP2=1;grade_PP2=3"
            continue

        # PP3 = 1
        if "PP3=1" in ev:

            if "grade_PP3=3" in ev:
                if "PP2=" not in ev:
                    evidence[v.key] += ";PP2=1;grade_PP2=3"

            elif "grade_PP3=2" in ev:
                if "PM1" not in labels:
                    if "PP2=" not in ev:
                        evidence[v.key] += ";PP2=1;grade_PP2=3"
                else:
                    if "PP2=" not in ev:
                        evidence[v.key] += ";PP2=0"

            elif "grade_PP3=1" in ev:
                if "PP2=" not in ev:
                    evidence[v.key] += ";PP2=0"

            else:
                if "PP2=" not in ev:
                    evidence[v.key] += ";PP2=0"

        else:
            if "PP2=" not in ev:
                evidence[v.key] += ";PP2=0"

    # ---------------- Script 10: BP7 ----------------
    for v in variants:
        if 'synonymous_variant' in v.Effect:
            if max(v.SpliceAI_acc_gain, v.SpliceAI_acc_loss, v.SpliceAI_don_gain, v.SpliceAI_don_loss) < 0.1:
                evidence[v.key] += ';BP7=1'

    # ---------------- Script 11: PP5 / BP6 ----------------
    for v in variants:
        evidence[v.key] += ';PP5=0;BP6=0'

    # ---------------- Script 12: PVS1 anula PP3 ----------------
    for v in variants:
        if 'PVS1' in v.InterVar_labels:
            evidence[v.key] += ';PP3=0'

    write_intervar(evidence, output_intervar_adjustments)

    return patients, maximum_credible_population_AF

# ============================================================
# CLI
# ============================================================

def parse_args():
    p = argparse.ArgumentParser(description="ACMG Bayesian annotation pipeline")
    p.add_argument("--variants", required=True, help="Excel with annotated variants")
    p.add_argument("--gene-inh", required=True, help="TSV with heritance models")
    p.add_argument("--gnomad", required=True, help="TSV constraint gnomAD")
    p.add_argument("--transcripts", required=True, help="JSON transcripts")
    p.add_argument("--output", required=True, help="File intervar_adjustments.tsv")
    p.add_argument("--prevalence", type=float, required=True)
    p.add_argument("--penetrance", type=float, required=True)
    p.add_argument("--inheritance_model", type=str, required=True)
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()
    run_pipeline(
        path_excel_variantes=args.variants,
        path_gene_inh_tsv=args.gene_inh,
        path_gnomad_constraint=args.gnomad,
        path_transcripts_json=args.transcripts,
        output_intervar_adjustments=args.output,
        prevalence=args.prevalence,
        penetrance=args.penetrance,
        inheritance_model=args.inheritance_model
    )

