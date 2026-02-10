# ACMG/ClinGen Variant Interpretation Pipeline for Congenital Hypogonadotropic Hypogonadism

This repository contains a reproducible bioinformatics pipeline designed to support
variant interpretation under ACMG/ClinGen guidelines in the context of congenital
hypogonadotropic hypogonadism (CHH).

The pipeline integrates curated patient-level variant data derived from a systematic
literature review with large-scale genomic annotations and Bayesian ACMG evidence
refinement.

This work was developed in the context of a systematic review performed at
Hospital de Niños Ricardo Gutiérrez.

---

## Pipeline Overview

The pipeline performs the following high-level operations:

1. Conversion of a curated variant database into VCF format
2. External variant annotation using standard resources
3. Bayesian ACMG scoring
4. Transcript-aware expansion of variant annotations
5. Integration with curated patient data using MANE transcripts
6. Iterative ACMG evidence refinement
7. Generation of a final refined variant database

## Scope and design considerations

This pipeline explicitly distinguishes between variant-level and transcript-level
representations:

- Transcript-level expansion is performed only when transcript-aware decisions are required
  (e.g. MANE selection, HGVS interpretation).
- Variant-level tables are used for Bayesian evidence aggregation and ACMG refinement steps.

This design avoids unintended duplication or loss of variants while maintaining
biological interpretability.

---

## Requirements

### Software

- Python ≥ 3.8
- Java ≥ 8 (for GATK tools)
- Access to external variant annotation pipelines

### External annotation resources (user-provided)

The annotated VCF must include at least:

- SnpEff
- gnomAD v4.1.0
- ClinVar
- REVEL
- SpliceAI
- InterVar

### Python dependencies

- pandas
- numpy
- openpyxl

---

## Input Files

| File | Description |
|-----|------------|
| `rs_base_18_3_24_raw.xlsx` | Curated patient and variant database |
| `transcripts.json` | MANE transcript definitions |
| `gene_inh_omim.tsv` | Gene inheritance modes |
| `gnomad.v4.1.constraint_metrics.tsv` | Gene-level constraint metrics |

---

## Step-by-Step Pipeline Execution

### Step 1 — Convert curated database to VCF

This step extracts variants from the curated database and converts them to VCF format.

```bash
python database_variants_2_vcf.py \
  --input syst_rev_test_raw.xlsx \
  --output variants.vcf \
  --sample SAMPLE
```

### Step 2 — External variant annotation

The generated VCF must be annotated using institutional or local bioinformatics
infrastructure.

Variant effect prediction is performed using SnpEff:

https://github.com/pcingola/SnpEff

Genomic coordinates and clinical significance are integrated from gnomAD and ClinVar using SnpSift. Functional and splicing impact scores are further derived from REVEL and SpliceAI via the dbNSFP plugin.

Variants are classified using InterVar:

https://github.com/WGLab/InterVar

The output of this step is an annotated VCF containing all required annotations.

### Step 3 — Initial Bayesian ACMG scoring

```bash
python bayes_2025.py \
  --vcf variants_annotated.vcf \
  --out variants_annotated_bayesian_initial.vcf
```

### Step 4 — Expand variants by transcript

Variants are expanded so that each transcript annotated by SnpEff (ANN field)
corresponds to a separate row.

```bash
python vcf_2_xlsx.py \
  --vcf variants_annotated_bayesian_initial.vcf \
  --out annotation_expanded_by_transcript.xlsx
```

### Step 5 — Merge expanded annotations with curated database

Annotations are merged back into the curated database using predefined MANE transcripts.

```bash
python merge_database_annotation.py \
  --base syst_rev_test_raw.xlsx \
  --annotation annotation_expanded_by_transcript.xlsx \
  --transcripts transcripts.json \
  --output bdd_annotated.xlsx \
  --annot-filtered annotation_expanded_by_transcript_filtered.xlsx
```

### Step 6 — ACMG evidence refinement (Step 1)

This step generates disease-specific ACMG evidence adjustments based on prevalence,
penetrance, inheritance model, and gene-level constraints.

```bash
python run_acmg_pipeline_step_1.py \
  --variants bdd_annotated.xlsx \
  --gene-inh gene_inh_omim.tsv \
  --gnomad gnomad.v4.1.constraint_metrics.tsv \
  --transcripts transcripts.json \
  --output intervar_adjustments_step_1.tsv \
  --prevalence 0.00025 \
  --penetrance 0.5 \
  --inheritance_model AR
```

### Step 7 — Bayesian update with refined evidence (Step 1)

```bash
python bayes_2025.py \
  --vcf variants_annotated.vcf \
  --out bayesian_step_1.vcf \
  --evidence intervar_adjustments_step_1.tsv
```

### Step 8 – Convert Bayesian-annotated VCF to tabular format (Step 1)

This step converts the Bayesian-adjusted VCF into a TSV table using GATK `VariantsToTable`.
The resulting table contains one row per variant and is used as input for subsequent ACMG
evidence refinement steps.

This step can be run locally on any system with Java and GATK installed.

```bash
gatk VariantsToTable   -V bayesian_step_1.vcf   -O table_bayesian_step_1.tsv   -F CHROM   -F POS   -F ID   -F REF   -F ALT   -F SNPEFF_ALLELE   -F SNPEFF_EXON_ID   -F SNPEFF_EXON_COUNT   -F SNPEFF_AMINO_ACID_CHANGE   -F SNPEFF_CODON_CHANGE   -F SNPEFF_EFFECT   -F SNPEFF_GENE_NAME   -F SNPEFF_IMPACT   -F SNPEFF_TRANSCRIPT_ID   -F CLNSIG   -F CLNREVSTAT   -F InterVarEvidence   -F InterVarVeredict   -F REVEL_score   -F gnomAD_exome_ALL   -F gnomAD_genome_ALL   -F gnomad312_AF_faf95_popmax   -F SAI_DS_AG   -F SAI_DS_AL   -F SAI_DS_DG   -F SAI_DS_DL   -F ACMG_BAYESIAN_SCORE   -F ACMG_BAYESIAN_VEREDICT   -F ANN
```

### Step 9 — ACMG evidence refinement (Step 2)

```bash
python run_acmg_pipeline_step_2.py \
  --excel bdd_annotated.xlsx \
  --intervar intervar_adjustments_step_1.tsv \
  --table-bayes table_bayesian_step_1.tsv \
  --output intervar_adjustments_step_2.tsv
```

### Step 10 — Bayesian update with refined evidence (Step 2)

```bash
python bayes_2025.py \
  --vcf variants_annotated.vcf \
  --out bayesian_step_2.vcf \
  --evidence intervar_adjustments_step_2.tsv
```

### Step 11 – Convert Bayesian-annotated VCF to tabular format (Step 2)

```bash
gatk VariantsToTable   -V bayesian_step_2.vcf   -O table_bayesian_step_2.tsv   -F CHROM   -F POS   -F ID   -F REF   -F ALT   -F SNPEFF_ALLELE   -F SNPEFF_EXON_ID   -F SNPEFF_EXON_COUNT   -F SNPEFF_AMINO_ACID_CHANGE   -F SNPEFF_CODON_CHANGE   -F SNPEFF_EFFECT   -F SNPEFF_GENE_NAME   -F SNPEFF_IMPACT   -F SNPEFF_TRANSCRIPT_ID   -F CLNSIG   -F CLNREVSTAT   -F InterVarEvidence   -F InterVarVeredict   -F REVEL_score   -F gnomAD_exome_ALL   -F gnomAD_genome_ALL   -F gnomad312_AF_faf95_popmax   -F SAI_DS_AG   -F SAI_DS_AL   -F SAI_DS_DG   -F SAI_DS_DL   -F ACMG_BAYESIAN_SCORE   -F ACMG_BAYESIAN_VEREDICT   -F ANN
```

### Step 12 — ACMG evidence refinement (Step 3)

MCPAF stands for Maximum Credible Population Allele Frequency. This value is displayed during execution of `run_acmg_pipeline_step_1.py`. In this example, the 0.000025 is given as an example.

```bash
python run_acmg_pipeline_step_3.py \
  --excel bdd_annotated.xlsx \
  --intervar intervar_adjustments_step_2.tsv \
  --table-bayes table_bayesian_step_2.tsv \
  --output intervar_adjustments_step_3.tsv \
  --gene-inh gene_inh_omim.tsv \
  --mcpaf 0.000025
```

### Step 13 — Bayesian update with refined evidence (Step 3)

```bash
python bayes_2025.py \
  --vcf variants_annotated.vcf \
  --out bayesian_step_3.vcf \
  --evidence intervar_adjustments_step_3.tsv
```

### Step 14 – Convert Bayesian-annotated VCF to tabular format (Step 3)

```bash
gatk VariantsToTable   -V bayesian_step_3.vcf   -O table_bayesian_step_3.tsv   -F CHROM   -F POS   -F ID   -F REF   -F ALT   -F SNPEFF_ALLELE   -F SNPEFF_EXON_ID   -F SNPEFF_EXON_COUNT   -F SNPEFF_AMINO_ACID_CHANGE   -F SNPEFF_CODON_CHANGE   -F SNPEFF_EFFECT   -F SNPEFF_GENE_NAME   -F SNPEFF_IMPACT   -F SNPEFF_TRANSCRIPT_ID   -F CLNSIG   -F CLNREVSTAT   -F InterVarEvidence   -F InterVarVeredict   -F REVEL_score   -F gnomAD_exome_ALL   -F gnomAD_genome_ALL   -F gnomad312_AF_faf95_popmax   -F SAI_DS_AG   -F SAI_DS_AL   -F SAI_DS_DG   -F SAI_DS_DL   -F ACMG_BAYESIAN_SCORE   -F ACMG_BAYESIAN_VEREDICT   -F ANN
```

### Step 15 — Final database refinement

```bash
python refine_database.py \
  --excel bdd_annotated.xlsx \
  --table-bayes table_bayesian_step_3.tsv \
  --output-tsv bdd_refined.tsv \
  --output-xlsx bdd_refined.xlsx
```

## Output

The final output is a refined variant database including:

- Curated patient information

- Transcript-aware variant annotations

- ACMG/ClinGen evidence assignments

- Bayesian pathogenicity classifications

## Notes

This pipeline is intended for research use only.

Clinical interpretation must be performed by qualified professionals.

Some steps require institutional computational resources.

No variant-level deduplication is applied prior to canonical transcript selection.
Each variant is guaranteed to be retained throughout the pipeline, ensuring
complete traceability from input to final output.

## Citation

Please cite the associated publication when using this pipeline.

DOI: 10.1093/hmg/ddag007

## Contact

For questions, issues, or suggestions, please use the GitHub issue tracker.
