#!/usr/bin/env python3
"""
Convert a curated Excel variant table into a single-sample VCF file.

This script is intended as Step 1 of the ACMG/ClinGen variant interpretation
pipeline. It reads a spreadsheet containing genomic variants and produces a
VCF v4.2 compliant file that can be further annotated by external tools.

Expected input columns (minimum):
- Chr : chromosome (e.g. 1, 2, X)
- Pos : 1-based genomic position
- Ref : reference allele
- Alt : alternate allele
- cigocity : zygosity (Hetero, Homo, Hemi)

Author: Franco Ginobbi
"""

import argparse
import pandas as pd
import vcfpy

# -----------------------------------------------------------------------------
# VCF header construction
# -----------------------------------------------------------------------------

def build_vcf_header(sample_name: str) -> vcfpy.Header:
    """Build a minimal but valid VCF header for a single-sample VCF."""

    samples = vcfpy.SamplesInfos([sample_name])

    header_lines = [
        vcfpy.header.HeaderLine("fileformat", "VCFv4.2"),
        vcfpy.header.FormatHeaderLine.from_mapping(
            vcfpy.OrderedDict([
                ("ID", "GT"),
                ("Number", 1),
                ("Type", "String"),
                ("Description", "Genotype"),
            ])
        ),
    ]

    # GRCh38 contig definitions
    contigs = {
        "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
        "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
        "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
        "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345,
        "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
        "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415,
    }

    for chrom, length in contigs.items():
        header_lines.append(
            vcfpy.header.ContigHeaderLine.from_mapping(
                vcfpy.OrderedDict([
                    ("ID", chrom),
                    ("length", length),
                    ("assembly", "GRCh38"),
                ])
            )
        )

    return vcfpy.Header(lines=header_lines, samples=samples)


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

def zygosity_to_gt(zygosity: str) -> str:
    """Convert zygosity label to VCF GT field."""

    if zygosity == "Hetero":
        return "0/1"
    if zygosity in {"Homo", "Hemi"}:
        return "1/1"
    return "./."


def infer_variant_type(ref: str, alt: str) -> str:
    """Infer a simple variant type label for ALT field."""

    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt) > 1:
        return "MNV"
    if len(ref) > len(alt) == 1:
        return "DEL"
    if len(ref) < len(alt) == 1:
        return "INS"
    return "INDEL"


# -----------------------------------------------------------------------------
# Core logic
# -----------------------------------------------------------------------------

def excel_to_vcf(input_excel: str, output_vcf: str, sample_name: str) -> None:
    """Convert Excel variant table to VCF."""

    df = pd.read_excel(input_excel)

    required_cols = {"Chr", "Pos", "Ref", "Alt", "cigocity"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    header = build_vcf_header(sample_name)
    writer = vcfpy.Writer.from_path(output_vcf, header)

    O = vcfpy.OrderedDict

    for _, row in df.iterrows():
        gt = zygosity_to_gt(str(row.cigocity))

        call = [vcfpy.Call(sample_name, O(GT=gt))]

        var_type = infer_variant_type(str(row.Ref), str(row.Alt))
        alt = vcfpy.Substitution(var_type, str(row.Alt))

        record = vcfpy.Record(
            CHROM=str(row.Chr),
            POS=int(row.Pos),
            ID=[],
            REF=str(row.Ref),
            ALT=[alt],
            QUAL=None,
            FILTER=[],
            INFO=O(),
            FORMAT=["GT"],
            calls=call,
        )

        writer.write_record(record)

    writer.close()


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert an Excel variant table into a single-sample VCF"
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input Excel file (.xlsx)",
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output VCF file",
    )

    parser.add_argument(
        "-s", "--sample",
        default="SAMPLE",
        help="Sample name (default: SAMPLE)",
    )

    args = parser.parse_args()

    excel_to_vcf(
        input_excel=args.input,
        output_vcf=args.output,
        sample_name=args.sample,
    )

    print(f"VCF successfully written to: {args.output}")


if __name__ == "__main__":
    main()

