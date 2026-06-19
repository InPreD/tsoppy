from contextlib import nullcontext
from os import path

from numpy import False_, float32, int32, str_
from polars import DataFrame
from pytest import mark, raises

from tsoppy.general.classes import (
    SmallVariantGenomeVcf,
    TmbTraceTsv,
    VariantsAnnotatedJson,
    WorkflowOutput,
)

# Define path to test data - cannot be absolute due to different paths locally and in CI
test_data_dir = "tests/test_data/general_classes"


@mark.parametrize(
    "inputs, exception, want",
    [
        (
            ("config.yaml", path.join(test_data_dir, "dragen/standard")),
            nullcontext(),
            "dragen_2.6.2.4",
        ),
        (
            ("config.yaml", path.join(test_data_dir, "localapp/standard")),
            nullcontext(),
            "localapp_ruo-2.2.0.12",
        ),
        (
            ("config.yaml", path.join(test_data_dir, "localapp/non-existent")),
            raises(FileNotFoundError),
            "",
        ),
    ],
)
def test_workflowoutput_init(inputs, exception, want):
    with exception:
        got = WorkflowOutput(inputs[0], inputs[1])
        assert got.workflow_id() == want


@mark.parametrize(
    "inputs, exception, want",
    [
        (
            ("config.yaml", path.join(test_data_dir, "dragen/standard"), "sample1"),
            nullcontext(),
            (
                "chr1",
                1000000,
                "A",
                int32(0),
                int32(700),
                int32(1),
                False_,
                float32(-1),
                str_("A/A"),
            ),
        ),
        (
            ("config.yaml", path.join(test_data_dir, "localapp/standard"), "sample1"),
            nullcontext(),
            (
                "chr1",
                1000000,
                "A",
                int32(0),
                int32(-1),
                int32(-1),
                False_,
                float32(0),
                str_("A/A"),
            ),
        ),
        (
            ("config.yaml", path.join(test_data_dir, "dragen/non-existent"), "sample1"),
            raises(FileNotFoundError),
            (
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            ),
        ),
    ],
)
def test_smallvariantgenomevcf_create(inputs, exception, want):
    with exception:
        workflow_output = WorkflowOutput(inputs[0], inputs[1])
        got = SmallVariantGenomeVcf.create(workflow_output, inputs[2])
        got_variants = list(got.vcf)
        assert len(got_variants) == 1
        got_variant = got_variants[0]
        assert got_variant.CHROM == want[0]
        assert got_variant.POS == want[1]
        assert got_variant.REF == want[2]
        assert got_variant.gt_types[0] == want[3]
        assert got_variant.gt_ref_depths[0] == want[4]
        assert got_variant.gt_alt_depths[0] == want[5]
        assert got_variant.gt_phases[0] == want[6]
        assert got_variant.gt_quals[0] == want[7]
        assert got_variant.gt_bases[0] == want[8]


@mark.parametrize(
    "inputs, exception, want",
    [
        (
            ("config.yaml", path.join(test_data_dir, "dragen/standard"), "sample1"),
            nullcontext(),
            DataFrame(
                {
                    "Chromosome": ["chr1"],
                    "Position": [1000000],
                    "RefCall": ["A"],
                    "AltCall": ["T"],
                    "VAF": [0.5],
                    "Depth": [550],
                    "CytoBand": ["1p1.1"],
                    "GeneName": ["GEN1"],
                    "VariantType": ["SNV"],
                    "CosmicIDs": ["COSM0001;COSM0002"],
                    "MaxCosmicCount": [2],
                    "ClinVarIDs": ["RCV0001.1"],
                    "ClinVarSignificance": ["not provided"],
                    "AlleleCountsGnomadExome": [100],
                    "AlleleCountsGnomadGenome": [100],
                    "AlleleCounts1000Genomes": [30],
                    "MaxDatabaseAlleleCounts": [100],
                    "GermlineFilterDatabase": [True],
                    "GermlineFilterProxi": [False],
                    "Nonsynonymous": [True],
                    "withinValidTmbRegion": [True],
                    "IncludedInTMBNumerator": [False],
                    "Status": ["Germline_DB"],
                    "ProteinChange": ["NP_001.2:p.(Ala1Ile)"],
                    "CDSChange": ["NM_000001.1:c.100A>T"],
                    "Exons": ["1/2"],
                    "Consequence": ["missense_variant"],
                }
            ),
        ),
        (
            ("config.yaml", path.join(test_data_dir, "localapp/standard"), "sample1"),
            nullcontext(),
            DataFrame(
                {
                    "Chromosome": ["chr1"],
                    "Position": [1000000],
                    "RefCall": ["A"],
                    "AltCall": ["T"],
                    "VAF": [0.5],
                    "Depth": [600],
                    "CytoBand": ["1p1.1"],
                    "GeneName": ["GEN1"],
                    "VariantType": ["SNV"],
                    "CosmicIDs": ["COSM0001;COSM0002"],
                    "MaxCosmicCount": [2],
                    "AlleleCountsGnomadExome": [100],
                    "AlleleCountsGnomadGenome": [100],
                    "AlleleCounts1000Genomes": [30],
                    "MaxDatabaseAlleleCounts": [100],
                    "GermlineFilterDatabase": [True],
                    "GermlineFilterProxi": [False],
                    "CodingVariant": [True],
                    "Nonsynonymous": [True],
                    "IncludedInTMBNumerator": [False],
                }
            ),
        ),
        (
            ("config.yaml", path.join(test_data_dir, "dragen/non-existent"), "sample1"),
            raises(FileNotFoundError),
            None,
        ),
    ],
)
def test_tmbtracetsv_create(inputs, exception, want):
    with exception:
        workflow_output = WorkflowOutput(inputs[0], inputs[1])
        got = TmbTraceTsv.create(workflow_output, inputs[2])
        assert got.table.equals(want)


@mark.parametrize(
    "inputs, exception, want",
    [
        (
            ("config.yaml", path.join(test_data_dir, "dragen/standard"), "sample1"),
            nullcontext(),
            {"id": "sample1"},
        ),
        (
            ("config.yaml", path.join(test_data_dir, "localapp/standard"), "sample1"),
            nullcontext(),
            {"id": "sample1"},
        ),
        (
            ("config.yaml", path.join(test_data_dir, "dragen/non-existent"), "sample1"),
            raises(FileNotFoundError),
            None,
        ),
    ],
)
def test_variantsannotatedjson_create(inputs, exception, want):
    with exception:
        workflow_output = WorkflowOutput(inputs[0], inputs[1])
        got = VariantsAnnotatedJson.create(workflow_output, inputs[2])
        assert got.data == want
