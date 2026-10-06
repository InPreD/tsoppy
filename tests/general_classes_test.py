from contextlib import nullcontext
from os import path

from numpy import False_, float32, int32, str_
from polars import DataFrame
from polars.testing import assert_frame_equal
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
        got = WorkflowOutput(*inputs)
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
                {
                    "fileformat": "VCFv4.2",
                    "FILTER": {"PASS": {"Description": '"All filters passed"'}},
                    "ALT": {
                        "NON_REF": {
                            "Description": '"Represents any possible alternative allele at this location"'
                        }
                    },
                    "FORMAT": {
                        "AD": {
                            "Number": "R",
                            "Type": "Integer",
                            "Description": '"Allelic depths (counting only informative reads out of the total reads) for the ref and alt alleles in the order listed"',
                        },
                        "DP": {
                            "Number": "1",
                            "Type": "Integer",
                            "Description": '"Approximate read depth (reads with MQ=255 or with bad mates are filtered)"',
                        },
                        "GT": {
                            "Number": "1",
                            "Type": "String",
                            "Description": '"Genotype"',
                        },
                        "SQ": {
                            "Number": "A",
                            "Type": "Float",
                            "Description": '"Somatic quality"',
                        },
                    },
                    "contig": {"chr1": {"length": "249250621"}},
                    "reference": "file://hashtable/reference.bin",
                },
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
                {
                    "fileformat": "VCFv4.1",
                    "FILTER": {"PASS": {"Description": '"All filters passed"'}},
                    "ALT": {
                        "NON_REF": {
                            "Description": '"Represents any possible alternative allele at this location"'
                        }
                    },
                    "FORMAT": {
                        "AD": {
                            "Number": None,
                            "Type": "Integer",
                            "Description": '"Allele Depth"',
                        },
                        "DP": {
                            "Number": "1",
                            "Type": "Integer",
                            "Description": '"Total Depth Used For Variant Calling"',
                        },
                        "GT": {
                            "Number": "1",
                            "Type": "String",
                            "Description": '"Genotype"',
                        },
                        "SQ": {
                            "Number": "A",
                            "Type": "Float",
                            "Description": '"Somatic quality"',
                        },
                        "GQ": {
                            "Number": "1",
                            "Type": "Integer",
                            "Description": '"Genotype Quality"',
                        },
                        "VF": {
                            "Number": None,
                            "Type": "Float",
                            "Description": '"Variant Frequency"',
                        },
                        "NL": {
                            "Number": "1",
                            "Type": "Integer",
                            "Description": '"Applied BaseCall Noise Level"',
                        },
                        "SB": {
                            "Number": "1",
                            "Type": "Float",
                            "Description": '"StrandBias Score"',
                        },
                        "NC": {
                            "Number": "1",
                            "Type": "Float",
                            "Description": '"Fraction of bases which were uncalled or with basecall quality below the minimum threshold"',
                        },
                        "US": {
                            "Number": None,
                            "Type": "Integer",
                            "Description": '"Supporting read type counts"',
                        },
                    },
                    "contig": {"chr1": {"length": "249250621"}},
                    "reference": "/opt/illumina/resources/genomes/hg19_hardPAR",
                },
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
                None,
            ),
        ),
    ],
)
def test_smallvariantgenomevcf_create(inputs, exception, want):
    with exception:
        workflow_output = WorkflowOutput(*inputs[:2])
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
        assert got.header_dict == want[9]


@mark.parametrize(
    "inputs, exception, want",
    [
        (
            ("config.yaml", path.join(test_data_dir, "dragen/standard"), "sample2"),
            nullcontext(),
            DataFrame(
                {
                    "variant_ID": [
                        "chr1:1000000:A>G",
                        "chr1:1000000:A>C",
                        "chr1:1200000:CTGAAAGAGC>GTGAAAGAGT",
                        "chr1:2000000:A>G",
                        "chr1:3000000:G>A",
                        "chr1:5000000:G>A",
                        "chr1:6000000:G>A",
                        "chr1:7000000:C>T",
                        "chr1:8000000:T>TC",
                        "chr1:9000000:CG>C",
                    ],
                    "variant_type": [
                        "SNV",
                        "SNV",
                        "MNV",
                        "SNV",
                        "SNV",
                        "SNV",
                        "SNV",
                        "SNV",
                        "insertion",
                        "deletion",
                    ],
                    "Illumina_variant_class": [
                        "VCF_filtered",
                        "Blacklist",
                        "Somatic",
                        "Germline_DB",
                        "VCF_filtered",
                        "Somatic",
                        "Germline_DB",
                        "Germline_Proxi",
                        "VCF_filtered",
                        "VCF_filtered",
                    ],
                    "overlapping_genes": [
                        ["SYN_GENE_A"],
                        ["GENE_4"],
                        ["SYN_GENE_X"],
                        ["SYN_GENE_A"],
                        ["SYN_GENE_B"],
                        ["SYN_GENE_B"],
                        ["SYN_GENE_C"],
                        ["SYN_GENE_D"],
                        ["SYN_GENE_E"],
                        ["SYN_GENE_F"],
                    ],
                    "AD": [
                        [189, 1, 0],
                        [177, 194, 0],
                        [244, 12, 0],
                        [419, 19, 0],
                        [0, 123, 0],
                        [0, 546, 0],
                        [0, 1243, 0],
                        [0, 186, 0],
                        [218, 26, 0],
                        [801, 6, 0],
                    ],
                    "DP": [190, 371, 256, 438, 123, 546, 1243, 186, 244, 807],
                    "AF": [
                        [0.0054, 0.0],
                        [0.5064, 0.0],
                        [0.042, 0.0],
                        [0.0434, 0.0],
                        [0.0150, 0.0],
                        [0.7738, 0.0],
                        [0.4788, 0.0],
                        [1.0, 0.0],
                        [0.0864, 0.0],
                        [0.0079, 0.0],
                    ],
                }
            ),
        ),
        (
            ("config.yaml", path.join(test_data_dir, "localapp/standard"), "sample2"),
            nullcontext(),
            DataFrame(
                {
                    "variant_ID": [
                        "chr1:1000000:A>G",
                        "chr1:1000000:A>C",
                        "chr1:1200000:CTGAAAGAGC>GTGAAAGAGT",
                        "chr1:2000000:A>G",
                        "chr1:3000000:G>A",
                        "chr1:5000000:G>A",
                        "chr1:6000000:G>A",
                        "chr1:7000000:C>T",
                        "chr1:8000000:T>TC",
                        "chr1:9000000:CG>C",
                    ],
                    "variant_type": [
                        "SNV",
                        "SNV",
                        "MNV",
                        "SNV",
                        "SNV",
                        "SNV",
                        "SNV",
                        "SNV",
                        "insertion",
                        "deletion",
                    ],
                    "Illumina_variant_class": [
                        "VCF_filtered",
                        "Blacklist",
                        "Somatic",
                        "Germline_DB",
                        "VCF_filtered",
                        "Somatic",
                        "Germline_DB",
                        "Germline_Proxi",
                        "VCF_filtered",
                        "VCF_filtered",
                    ],
                    "overlapping_genes": [
                        ["SYN_GENE_A"],
                        ["GENE_4"],
                        ["SYN_GENE_X"],
                        ["SYN_GENE_A"],
                        ["SYN_GENE_B"],
                        ["SYN_GENE_B"],
                        ["SYN_GENE_C"],
                        ["SYN_GENE_D"],
                        ["SYN_GENE_E"],
                        ["SYN_GENE_F"],
                    ],
                    "AD": [
                        [189, 1],
                        [177, 194],
                        [244, 12],
                        [419, 19],
                        [0, 123],
                        [0, 546],
                        [0, 1243],
                        [0, 186],
                        [218, 26],
                        [801, 6],
                    ],
                    "DP": [190, 371, 256, 438, 123, 546, 1243, 186, 244, 807],
                    "VF": [
                        0.0054,
                        0.5064,
                        0.042,
                        0.0434,
                        0.0150,
                        0.7738,
                        0.4788,
                        1.0,
                        0.0864,
                        0.0079,
                    ],
                }
            ),
        ),
    ],
)
def test_smallvariantgenomevcf_merge(inputs, exception, want):
    with exception:
        workflow_output = WorkflowOutput(*inputs[:2])
        vcf_obj = SmallVariantGenomeVcf.create(workflow_output, inputs[2])
        _ = vcf_obj.merge(workflow_output, inputs[2])
        assert_frame_equal(left=vcf_obj.merged_table, right=want, check_dtypes=False)


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
        workflow_output = WorkflowOutput(*inputs[:2])
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
        workflow_output = WorkflowOutput(*inputs[:2])
        got = VariantsAnnotatedJson.create(workflow_output, inputs[2])
        assert got.data == want
