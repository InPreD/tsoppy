"""Unit tests for the add_variant_summary main module."""

import polars
from pytest import mark

from tsoppy.add_variant_summary.main import (
    _detect_data_type,
    _format_cnv,
    _format_dna_row,
    _format_fusions,
    _format_gene_amp,
    _format_large_rearr,
    _format_loh,
    _format_msi,
    _format_rna_row,
    _format_splice_vars,
    _format_tmb,
    _parse_sample,
)


# --- Helpers to build transposed key-value and tabular DataFrames ----------


def _analysis_details(dna="IPD1111-D", rna="NA"):
    return polars.DataFrame({"DNA Sample ID": [dna], "RNA Sample ID": [rna]})


def _tmb(total="7.8", eligible="10"):
    return polars.DataFrame(
        {"Total TMB": [total], "Number of Passing Eligible Variants": [eligible]}
    )


def _msi(percent="4.13", unstable="5", usable="121"):
    return polars.DataFrame(
        {
            "Percent Unstable MSI Sites": [percent],
            "Total MSI Sites Unstable": [unstable],
            "Usable MSI Sites": [usable],
        }
    )


def _gis(score="32", tumor_fraction="0.650", ploidy="3.92"):
    return polars.DataFrame(
        {
            "Genomic Instability Score": [score],
            "Tumor Fraction": [tumor_fraction],
            "Ploidy": [ploidy],
        }
    )


def _cnv():
    return polars.DataFrame(
        {
            "Gene": ["EGFR", "MET"],
            "Fold Change": ["1.5", "2.0"],
            "Copy Number Variant": ["DUP", "DEL"],
            "Absolute Copy Number": ["5", "2"],
        }
    )


def _loh():
    return polars.DataFrame(
        {"Gene": ["BRCA1"], "Absolute Copy Number": ["1"], "Minor Copy Number": ["0"]}
    )


def _gene_amp():
    return polars.DataFrame({"Gene": ["MDM4"], "Fold Change": ["1.741"]})


def _splice():
    return polars.DataFrame(
        {
            "Gene": ["EGFR"],
            "Affected Exon": ["2-3/28"],
            "Breakpoint 1": ["chr7:55087058"],
            "Breakpoint 2": ["chr7:55214298"],
            "Splice Supporting Reads": ["15"],
            "Reference Reads Transcript": ["13089"],
        }
    )


def _fusions():
    return polars.DataFrame(
        {
            "Gene Pair": ["USP13-PIK3CA"],
            "Breakpoint 1": ["chr3:179463005"],
            "Breakpoint 2": ["chr3:178947060"],
            "Fusion Supporting Reads": ["18"],
            "Gene 1 Reference Reads": ["0"],
            "Gene 2 Reference Reads": ["429"],
        }
    )


def _large_rearr():
    return polars.DataFrame({"Gene": ["BRCA1", "BRCA2"]})


def _dragen_sections(dna="IPD1111-D", rna="NA"):
    return {
        "Analysis Details": _analysis_details(dna, rna),
        "TMB": _tmb(),
        "MSI": _msi(),
        "GIS": _gis(),
        "Copy Number Variants": _cnv(),
        "Loss of Heterozygosity": _loh(),
        "Splice Variants": _splice(),
        "Fusions": _fusions(),
        "Large Rearrangements": _large_rearr(),
    }


def _localapp_sections(dna="IPD1111-D", rna="NA"):
    return {
        "Analysis Details": _analysis_details(dna, rna),
        "TMB": _tmb(),
        "MSI": _msi(),
        "Gene Amplifications": _gene_amp(),
        "Splice Variants": _splice(),
        "Fusions": _fusions(),
    }


# --- Key-value formatters ----------------------------------------------------


@mark.parametrize(
    "sections, want",
    [
        ({"TMB": _tmb()}, "7.8 (10)"),
        ({"TMB": _tmb("NA", "10")}, "NA"),
        ({}, "NA"),
    ],
)
def test_format_tmb(sections, want):
    assert _format_tmb(sections) == want


@mark.parametrize(
    "sections, want",
    [
        ({"MSI": _msi()}, "4.13 (5/121)"),
        ({"MSI": _msi("NA", "5", "121")}, "NA"),
        ({}, "NA"),
    ],
)
def test_format_msi(sections, want):
    assert _format_msi(sections) == want


# --- Tabular formatters ------------------------------------------------------


@mark.parametrize(
    "sections, want",
    [
        ({"Copy Number Variants": _cnv()}, "EGFR (1.5, DUP, 5)|MET (2.0, DEL, 2)"),
        ({"Copy Number Variants": polars.DataFrame()}, "NA"),
        ({}, "NA"),
    ],
)
def test_format_cnv(sections, want):
    assert _format_cnv(sections) == want


@mark.parametrize(
    "sections, want",
    [
        ({"Loss of Heterozygosity": _loh()}, "BRCA1 (1, 0)"),
        ({"Loss of Heterozygosity": polars.DataFrame()}, "NA"),
        ({}, "NA"),
    ],
)
def test_format_loh(sections, want):
    assert _format_loh(sections) == want


@mark.parametrize(
    "sections, want",
    [
        ({"Gene Amplifications": _gene_amp()}, "MDM4 (1.741, <DUP>)"),
        ({"Gene Amplifications": polars.DataFrame()}, "NA"),
        ({}, "NA"),
    ],
)
def test_format_gene_amp(sections, want):
    assert _format_gene_amp(sections) == want


@mark.parametrize(
    "sections, want",
    [
        (
            {"Splice Variants": _splice()},
            "EGFR[2-3/28] (chr7:55087058-chr7:55214298) 15/13089",
        ),
        ({"Splice Variants": polars.DataFrame()}, "NA"),
        ({}, "NA"),
    ],
)
def test_format_splice_vars(sections, want):
    assert _format_splice_vars(sections) == want


@mark.parametrize(
    "sections, want",
    [
        (
            {"Fusions": _fusions()},
            "USP13-PIK3CA (chr3:179463005-chr3:178947060) 18/0/429",
        ),
        ({"Fusions": polars.DataFrame()}, "NA"),
        ({}, "NA"),
    ],
)
def test_format_fusions(sections, want):
    assert _format_fusions(sections) == want


@mark.parametrize(
    "sections, want",
    [
        ({"Large Rearrangements": _large_rearr()}, "BRCA1|BRCA2"),
        ({"Large Rearrangements": polars.DataFrame()}, "NA"),
        ({}, "NA"),
    ],
)
def test_format_large_rearr(sections, want):
    assert _format_large_rearr(sections) == want


# --- Data type detection -----------------------------------------------------


@mark.parametrize(
    "sections, want",
    [
        (_dragen_sections("IPD1111-D", "NA"), (True, False)),
        (_dragen_sections("NA", "IPD1111-R"), (False, True)),
        (_dragen_sections("IPD1111-D", "IPD1111-R"), (True, True)),
        (_dragen_sections("NA", "NA"), (False, False)),
    ],
)
def test_detect_data_type(sections, want):
    assert _detect_data_type(sections) == want


# --- Row builders ------------------------------------------------------------


@mark.parametrize(
    "sections, workflow_type, want",
    [
        (
            _dragen_sections(),
            "dragen",
            [
                "IPD1111-D",
                "7.8 (10)",
                "4.13 (5/121)",
                "EGFR (1.5, DUP, 5)|MET (2.0, DEL, 2)",
                "BRCA1 (1, 0)",
                "NA",  # splice_variants
                "NA",  # fusions
                "32",  # GIS
                "0.650",  # Tumor fraction
                "3.92",  # Ploidy
                "BRCA1|BRCA2",  # Large rearrangements
            ],
        ),
        (
            _localapp_sections(),
            "localapp",
            [
                "IPD1111-D",
                "7.8 (10)",
                "4.13 (5/121)",
                "MDM4 (1.741, <DUP>)",
                "NA",
                "NA",
                "NA",
                "NA",
                "NA",
                "NA",
                "NA",
            ],
        ),
    ],
)
def test_format_dna_row(sections, workflow_type, want):
    assert _format_dna_row(sections, workflow_type) == want


@mark.parametrize(
    "sections, want",
    [
        (
            _dragen_sections("NA", "IPD1111-R"),
            [
                "IPD1111-R",
                "NA",
                "NA",
                "NA",
                "NA",
                "EGFR[2-3/28] (chr7:55087058-chr7:55214298) 15/13089",
                "USP13-PIK3CA (chr3:179463005-chr3:178947060) 18/0/429",
                "NA",
                "NA",
                "NA",
                "NA",
            ],
        ),
    ],
)
def test_format_rna_row(sections, want):
    assert _format_rna_row(sections) == want


# --- Integration: one parsed file -> rows -----------------------------------


@mark.parametrize(
    "sections, workflow_type, want",
    [
        (
            _dragen_sections("IPD1111-D", "NA"),
            "dragen",
            [
                [
                    "IPD1111-D",
                    "7.8 (10)",
                    "4.13 (5/121)",
                    "EGFR (1.5, DUP, 5)|MET (2.0, DEL, 2)",
                    "BRCA1 (1, 0)",
                    "NA",
                    "NA",
                    "32",
                    "0.650",
                    "3.92",
                    "BRCA1|BRCA2",
                ],
            ],
        ),
        (
            _dragen_sections("IPD1111-D", "IPD1111-R"),
            "dragen",
            [
                [
                    "IPD1111-D",
                    "7.8 (10)",
                    "4.13 (5/121)",
                    "EGFR (1.5, DUP, 5)|MET (2.0, DEL, 2)",
                    "BRCA1 (1, 0)",
                    "NA",
                    "NA",
                    "32",
                    "0.650",
                    "3.92",
                    "BRCA1|BRCA2",
                ],
                [
                    "IPD1111-R",
                    "NA",
                    "NA",
                    "NA",
                    "NA",
                    "EGFR[2-3/28] (chr7:55087058-chr7:55214298) 15/13089",
                    "USP13-PIK3CA (chr3:179463005-chr3:178947060) 18/0/429",
                    "NA",
                    "NA",
                    "NA",
                    "NA",
                ],
            ],
        ),
        (
            _localapp_sections("IPD1111-D", "NA"),
            "localapp",
            [
                [
                    "IPD1111-D",
                    "7.8 (10)",
                    "4.13 (5/121)",
                    "MDM4 (1.741, <DUP>)",
                    "NA",
                    "NA",
                    "NA",
                    "NA",
                    "NA",
                    "NA",
                    "NA",
                ],
            ],
        ),
    ],
)
def test_parse_sample(sections, workflow_type, want):
    assert _parse_sample(sections, workflow_type) == want
