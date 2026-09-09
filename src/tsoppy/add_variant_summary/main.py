"""Create the variant summary table using the combined variant output files (cvo) for each sample."""

import logging
from pathlib import Path

import polars

from tsoppy.general.classes import CombinedVariantOutputTsv

logger = logging.getLogger(__name__)


# Set up variable names for data we want to extract
TMB_TOTAL = "Total TMB"
TMB_ELIGIBLE = "Number of Passing Eligible Variants"
MSI_PERCENT = "Percent Unstable MSI Sites"
MSI_UNSTABLE = "Total MSI Sites Unstable"
MSI_USABLE = "Usable MSI Sites"
GIS_SCORE = "Genomic Instability Score"
TUMOR_FRACTION = "Tumor Fraction"
PLOIDY = "Ploidy"

# Output table columns, in order.
HEADER = [
    "Sample_id",
    "TMB",
    "MSI",
    # dragen: Copy Number Variants; localapp: Gene Amplifications
    "Gene_amplifications_CNV",
    "LOH",
    "Splice_variants",
    "Fusions",
    "GIS",
    "Tumor_fraction",
    "Ploidy",
    "Large_rearrangements",
]


def Summarise_variant_output(
    config_yaml: str | Path,
    inpred_nomenclature: str | Path,
    results_directory: str | Path,
    output_file: str | Path,
) -> None:
    """Find all *_CombinedVariantOutput.tsv files in a directory tree and parse their results into a summary table."""
    cvo = CombinedVariantOutputTsv(config_yaml, inpred_nomenclature, results_directory)

    workflow_type = cvo.workflow_type  # detect workflow type
    with open(output_file, "w") as out_file:
        _write_header(out_file, cvo.headers)
        for file_path, sections in zip(cvo.files, cvo.sections):
            for row in _parse_sample(sections, workflow_type):
                out_file.write("\t".join(row) + "\n")


def _detect_data_type(sections: dict[str, polars.DataFrame]) -> tuple[bool, bool]:
    """Check the [Analysis Details] section DNA and RNA sample IDs to determine the file type."""
    details = sections.get("Analysis Details")
    has_dna = _value(details, "DNA Sample ID") != "NA"
    has_rna = _value(details, "RNA Sample ID") != "NA"
    return has_dna, has_rna


def _format_dna_row(
    sections: dict[str, polars.DataFrame], workflow_type: str
) -> list[str]:
    sample_id = _value(sections.get("Analysis Details"), "DNA Sample ID")

    if workflow_type == "dragen":
        gene_amp = _format_cnv(sections)
        loh = _format_loh(sections)
        gis = _value(sections.get("GIS"), GIS_SCORE)
        tumor_fraction = _value(sections.get("GIS"), TUMOR_FRACTION)
        ploidy = _value(sections.get("GIS"), PLOIDY)
        large_rearr = _format_large_rearr(sections)
    else:
        gene_amp = _format_gene_amp(sections)
        loh = "NA"
        gis = "NA"
        tumor_fraction = "NA"
        ploidy = "NA"
        large_rearr = "NA"

    return [
        sample_id,
        _format_tmb(sections),
        _format_msi(sections),
        gene_amp,
        loh,
        "NA",  # splice_variants (RNA only)
        "NA",  # fusions (RNA only)
        gis,
        tumor_fraction,
        ploidy,
        large_rearr,
    ]


def _format_rna_row(sections: dict[str, polars.DataFrame]) -> list[str]:
    sample_id = _value(sections.get("Analysis Details"), "RNA Sample ID")

    return [
        sample_id,
        "NA",  # TMB
        "NA",  # MSI
        "NA",  # gene_amplifications (DNA only)
        "NA",  # LOH (DNA only)
        _format_splice_vars(sections),
        _format_fusions(sections),
        "NA",  # GIS
        "NA",  # Tumor fraction
        "NA",  # Ploidy
        "NA",  # Large Rearrangements
    ]


def _write_header(out_file, headers: list[str]) -> None:
    for h in headers:
        out_file.write(f"# {h}\n")
    out_file.write("# Field formats:\n")
    out_file.write(
        "# - Tumor mutation burden (TMB): Total TMB (number of passing eligible variants)\n"
    )
    out_file.write(
        "# - Microsatellite instability (MSI): Percent Unstable MSI Sites (unstable/usable)\n"
    )
    out_file.write(
        "# - Gene amplifications (LocalApp) / Copy Number Variants (CNV) (DRAGEN V2): gene (fold_change, copy_number_variant, absolute_copy_number (DRAGEN V2 only))\n"
    )
    out_file.write(
        "# - Loss of Heterozygosity (LOH): gene (absolute_copy_number, minor_copy_number)\n"
    )
    out_file.write(
        "# - Splice variants: gene[affected_exon] (breakpoint_1-breakpoint_2) supporting/ref_reads\n"
    )
    out_file.write(
        "# - Fusions: gene_1-gene_2 (breakpoint_1-breakpoint_2) supporting/ref_gene_1/ref_gene_2\n"
    )
    out_file.write(
        "# - Genomic Instability Score (GIS) (DRAGEN V2 only): Genomic Instability Score, "
        "from the [GIS] section\n"
    )
    out_file.write(
        "# - Tumor fraction (DRAGEN V2 only): Tumor Fraction, from the [GIS] section\n"
    )
    out_file.write("# - Ploidy (DRAGEN V2 only): Ploidy, from the [GIS] section\n")
    out_file.write("# - Large Rearrangements: gene (DRAGEN V2 only)\n")
    out_file.write("\t".join(HEADER) + "\n")


def _parse_sample(
    sections: dict[str, polars.DataFrame], workflow_type: str
) -> list[list[str]]:
    has_dna, has_rna = _detect_data_type(sections)
    rows = []
    if has_dna:
        rows.append(_format_dna_row(sections, workflow_type))
    if has_rna:
        rows.append(_format_rna_row(sections))
    return rows


def _format_tmb(sections: dict[str, polars.DataFrame]) -> str:
    """Format TMB as "Total TMB (eligible_variants)"."""
    total_tmb = _value(sections.get("TMB"), TMB_TOTAL)
    if total_tmb == "NA":
        return "NA"
    return f"{total_tmb} ({_value(sections.get('TMB'), TMB_ELIGIBLE)})"


def _format_msi(sections: dict[str, polars.DataFrame]) -> str:
    """Format MSI as "percent (unstable/usable)"."""
    percent = _value(sections.get("MSI"), MSI_PERCENT)
    if percent == "NA":
        return "NA"
    return (
        f"{percent} ({_value(sections.get('MSI'), MSI_UNSTABLE)}/"
        f"{_value(sections.get('MSI'), MSI_USABLE)})"
    )


def _format_cnv(sections: dict[str, polars.DataFrame]) -> str:
    """Format dragen copy number variants from the [Copy Number Variants] section."""
    cna = sections.get("Copy Number Variants")
    if cna is None or cna.is_empty():
        return "NA"
    output = []
    for row in cna.iter_rows(named=True):
        if row["Gene"] in (None, "NA", "-") or "BETA FEATURE" in str(row["Gene"]):
            continue
        output.append(
            f"{row['Gene']} ({row['Fold Change']}, "
            f"{row['Copy Number Variant']}, {row['Absolute Copy Number']})"
        )
    return "|".join(output) if output else "NA"


def _format_loh(sections: dict[str, polars.DataFrame]) -> str:
    """Format loss of heterozygosity from the [Loss of Heterozygosity] section."""
    loh = sections.get("Loss of Heterozygosity")
    if loh is None or loh.is_empty():
        return "NA"
    output = []
    for row in loh.iter_rows(named=True):
        if row["Gene"] in (None, "NA", "-"):
            continue
        output.append(
            f"{row['Gene']} ({row['Absolute Copy Number']}, {row['Minor Copy Number']})"
        )
    return "|".join(output) if output else "NA"


def _format_splice_vars(sections: dict[str, polars.DataFrame]) -> str:
    """Format splice variants from the [Splice Variants] section."""
    splice = sections.get("Splice Variants")
    if splice is None or splice.is_empty():
        return "NA"
    output = []
    for row in splice.iter_rows(named=True):
        if row["Gene"] in (None, "NA", "-"):
            continue
        output.append(
            f"{row['Gene']}[{row['Affected Exon']}] "
            f"({row['Breakpoint 1']}-{row['Breakpoint 2']}) "
            f"{row['Splice Supporting Reads']}/{row['Reference Reads Transcript']}"
        )
    return "|".join(output) if output else "NA"


def _format_fusions(sections: dict[str, polars.DataFrame]) -> str:
    """Format fusions from the [Fusions] section."""
    fusions = sections.get("Fusions")
    if fusions is None or fusions.is_empty():
        return "NA"
    output = []
    for row in fusions.iter_rows(named=True):
        if row["Gene Pair"] in (None, "NA", "-"):
            continue
        output.append(
            f"{row['Gene Pair']} ({row['Breakpoint 1']}-{row['Breakpoint 2']}) "
            f"{row['Fusion Supporting Reads']}/{row['Gene 1 Reference Reads']}/"
            f"{row['Gene 2 Reference Reads']}"
        )
    return "|".join(output) if output else "NA"


def _format_large_rearr(sections: dict[str, polars.DataFrame]) -> str:
    """Format large rearrangements from the [Large Rearrangements] section."""
    lr = sections.get("Large Rearrangements")
    if lr is None or lr.is_empty():
        return "NA"
    output = []
    for row in lr.iter_rows(named=True):
        if row["Gene"] in (None, "NA", "-"):
            continue
        output.append(str(row["Gene"]))
    return "|".join(output) if output else "NA"


def _format_gene_amp(sections: dict[str, polars.DataFrame]) -> str:
    """Format localapp gene amplifications (all are <DUP>)."""
    gamp = sections.get("Gene Amplifications")
    if gamp is None or gamp.is_empty():
        return "NA"
    output = []
    for row in gamp.iter_rows(named=True):
        if row["Gene"] in (None, "NA", "-"):
            continue
        output.append(f"{row['Gene']} ({row['Fold Change']}, <DUP>)")
    return "|".join(output) if output else "NA"


def _value(df: polars.DataFrame | None, column: str, default: str = "NA") -> str:
    """Return the first value of a column, or `default` if missing/empty."""
    if df is None or df.is_empty() or column not in df.columns:
        return default
    value = df[column][0]
    return default if value is None else str(value)
