"""
update small variant vcf list subpackage main module unit tests.
"""

import filecmp
import os
import unittest
from contextlib import nullcontext
from os import path
from pytest import mark, raises

import polars

from tsoppy.update_small_variant_vcf_list.main import (
    VariantRecurrenceTable,
    InvalidSampleType,
    Vcf,
    VcfList,
)

# Define path to test data - cannot be absolute due to different paths locally and in CI
test_data_dir = "tests/test_data/update_small_variant_vcf_list_main"

# test constants
glob_pattern = "**/Results/**/*_MergedSmallVariants.genome.vcf"
tumor_sample_types = "C,D,d,L,M,P,p,R,r,T,X"
inpred_id_regex = (
    r"(?P<patient_id>\D{3}\d{4})-\D\d{2}-(?P<sample_type>\D)\d{2}-\D\d{2}.*.vcf$"
)

mock_variantrecurrencetable_info_header = [
    '# The recurrence summary values are written in format "X:A+B+C+D=M/N"',
    '#   - X: sample type, one of "T" = tumor, "N" = normal, "A" = any',
    "#   - A: number of samples of type X, in which given variant was seen with VAF < 0.01",
    "#   - B: number of samples of type X, in which given variant was seen with 0.01 <= VAF < 0.05",
    "#   - C: number of samples of type X, in which given variant was seen with 0.05 <= VAF < 0.35",
    "#   - D: number of samples of type X, in which given variant was seen with 0.35 <= VAF",
    "#   - M: number of samples of type X, in which given variant was seen with any VAF",
    "#   - N: number of investigated samples of type X, in which given variant was callable (i.e., the variant site had coverage >= 20)",
]


@mark.parametrize(
    "input, exception, want",
    [
        (
            path.join(test_data_dir, "variantrecurrencetable_init/non_existent.tsv"),
            nullcontext(),
            (
                polars.DataFrame({"sample_vcf": [], "sample_type": []}),
                mock_variantrecurrencetable_info_header,
                polars.DataFrame(
                    {
                        "variant_id": [],
                        "tumor_recurrence_summary": [],
                        "normal_recurrence_summary": [],
                        "total_recurrence_summary": [],
                    }
                ),
                False,
            ),
        ),
        (
            path.join(test_data_dir, "variantrecurrencetable_init/empty_file.tsv"),
            nullcontext(),
            (
                polars.DataFrame({"sample_vcf": [], "sample_type": []}),
                mock_variantrecurrencetable_info_header,
                polars.DataFrame(
                    {
                        "variant_id": [],
                        "tumor_recurrence_summary": [],
                        "normal_recurrence_summary": [],
                        "total_recurrence_summary": [],
                    }
                ),
                False,
            ),
        ),
        (
            path.join(
                test_data_dir, "variantrecurrencetable_init/successfully_parse_file.tsv"
            ),
            nullcontext(),
            (
                polars.DataFrame(
                    {
                        "sample_vcf": [
                            "/data/Logs_Intermediates/DnaDragenCaller/sample1/sample1_MergedSmallVariants.genome.vcf",
                            "/data/Logs_Intermediates/DnaDragenCaller/sample1/sample1_MergedSmallVariants.genome.vcf",
                        ],
                        "sample_type": ["N", "T"],
                    }
                ),
                mock_variantrecurrencetable_info_header,
                polars.DataFrame(
                    {
                        "variant_id": ["1:100:A>C"],
                        "tumor_recurrence_summary": ["T:1+0+0+0=1/1"],
                        "normal_recurrence_summary": ["N:0+0+0+0=0/0"],
                        "total_recurrence_summary": ["A:1+0+0+0=1/1"],
                    }
                ),
                True,
            ),
        ),
        ("non_existent/non_existent.tsv", raises(FileNotFoundError), None),
        (
            path.join(test_data_dir, "variantrecurrencetable_init/missing_column.tsv"),
            raises(ValueError),
            None,
        ),
        (
            path.join(
                test_data_dir, "variantrecurrencetable_init/duplicated_variant.tsv"
            ),
            raises(ValueError),
            None,
        ),
    ],
)
def test_variantrecurrencetable_init(input, exception, want):
    with exception:
        got = VariantRecurrenceTable(input)
        assert got.sample_vcf_header.equals(want[0])
        assert got.info_header == want[1]
        assert got.body.equals(want[2])
        assert got.update == want[3]


class TestVcfList(unittest.TestCase):
    def test_update(self):
        test_cases = [
            {
                "name": "successfully update small variant vcf list",
                "results_dir": path.join(
                    test_data_dir, "successfully_update_small_variant_vcf_list"
                ),
                "glob_pattern": glob_pattern,
                "vcf_list": path.join(
                    test_data_dir,
                    "successfully_update_small_variant_vcf_list/TSO500_vcf_list.tsv",
                ),
                "inpred_id_regex": inpred_id_regex,
                "tumor_sample_types": tumor_sample_types,
                "output": path.join(
                    test_data_dir,
                    "successfully_update_small_variant_vcf_list/TSO500_vcf_list_updated.tsv",
                ),
                "expected": path.join(
                    test_data_dir,
                    "successfully_update_small_variant_vcf_list/TSO500_vcf_list_expected.tsv",
                ),
            },
            {
                "name": "create new small variant vcf list",
                "results_dir": path.join(
                    test_data_dir, "create_new_small_variant_vcf_list"
                ),
                "glob_pattern": glob_pattern,
                "vcf_list": None,
                "inpred_id_regex": inpred_id_regex,
                "tumor_sample_types": tumor_sample_types,
                "output": path.join(
                    test_data_dir,
                    "create_new_small_variant_vcf_list/TSO500_vcf_list_updated.tsv",
                ),
                "expected": path.join(
                    test_data_dir,
                    "create_new_small_variant_vcf_list/TSO500_vcf_list_expected.tsv",
                ),
            },
            {
                "name": "small variant vcf list does not exist",
                "results_dir": path.join(
                    test_data_dir, "small_variant_vcf_list_does_not_exist"
                ),
                "glob_pattern": glob_pattern,
                "vcf_list": path.join(
                    test_data_dir,
                    "small_variant_vcf_list_does_not_exist/TSO500_vcf_list.tsv",
                ),
                "inpred_id_regex": inpred_id_regex,
                "tumor_sample_types": tumor_sample_types,
                "output": path.join(
                    test_data_dir,
                    "small_variant_vcf_list_does_not_exist/TSO500_vcf_list_updated.tsv",
                ),
                "expected": path.join(
                    test_data_dir,
                    "small_variant_vcf_list_does_not_exist/TSO500_vcf_list_expected.tsv",
                ),
            },
            {
                "name": "skip existing vcf",
                "results_dir": path.join(test_data_dir, "skip_existing_vcf"),
                "glob_pattern": glob_pattern,
                "vcf_list": path.join(
                    test_data_dir, "skip_existing_vcf/TSO500_vcf_list.tsv"
                ),
                "inpred_id_regex": inpred_id_regex,
                "tumor_sample_types": tumor_sample_types,
                "output": path.join(
                    test_data_dir, "skip_existing_vcf/TSO500_vcf_list_updated.tsv"
                ),
                "expected": path.join(
                    test_data_dir, "skip_existing_vcf/TSO500_vcf_list_expected.tsv"
                ),
            },
            {
                "name": "sample is control",
                "results_dir": path.join(test_data_dir, "sample_is_control"),
                "glob_pattern": glob_pattern,
                "vcf_list": path.join(
                    test_data_dir, "sample_is_control/TSO500_vcf_list.tsv"
                ),
                "inpred_id_regex": inpred_id_regex,
                "tumor_sample_types": tumor_sample_types,
                "output": path.join(
                    test_data_dir, "sample_is_control/TSO500_vcf_list_updated.tsv"
                ),
                "expected": path.join(
                    test_data_dir, "sample_is_control/TSO500_vcf_list_expected.tsv"
                ),
            },
        ]

        for test_case in test_cases:
            with self.subTest(msg=test_case["name"]):
                got = VcfList(
                    test_case["results_dir"],
                    test_case["glob_pattern"],
                    test_case["vcf_list"],
                    test_case["inpred_id_regex"],
                    test_case["tumor_sample_types"],
                    test_case["output"],
                )
                got.update()
                assert filecmp.cmp(test_case["output"], test_case["expected"])
                os.remove(test_case["output"])


class TestVcf(unittest.TestCase):
    def test_init(self):
        test_cases = [
            {
                "name": "include sample",
                "vcf": "IPH0001-D01-T01-A01_MergedSmallVariants.genome.vcf",
                "inpred_id_regex": inpred_id_regex,
                "tumor_sample_types": tumor_sample_types,
                "exception": nullcontext(),
                "patient_id": "IPH0001",
                "sample_type": "T",
            },
            {
                "name": "inpred id is not parsable",
                "vcf": "IPH0001D01-T01-A01_MergedSmallVariants.genome.vcf",
                "inpred_id_regex": inpred_id_regex,
                "tumor_sample_types": tumor_sample_types,
                "exception": raises(AttributeError),
                "patient_id": None,
                "sample_type": None,
            },
            {
                "name": "sample is neither tumor nor normal",
                "vcf": "IPH0001-D01-A01-A01_MergedSmallVariants.genome.vcf",
                "inpred_id_regex": inpred_id_regex,
                "tumor_sample_types": tumor_sample_types,
                "exception": raises(InvalidSampleType),
                "patient_id": "IPH0001",
                "sample_type": "A",
            },
        ]

        for test_case in test_cases:
            with self.subTest(msg=test_case["name"]):
                with test_case["exception"]:
                    got = Vcf(
                        test_case["vcf"],
                        test_case["inpred_id_regex"],
                        test_case["tumor_sample_types"],
                    )
                    assert got.patient_id == test_case["patient_id"]
                    assert got.sample_type == test_case["sample_type"]

    def test_row(self):
        test_cases = [
            {
                "name": "successfully return row",
                "vcf": "IPH0001-D01-T01-A01_MergedSmallVariants.genome.vcf",
                "inpred_id_regex": inpred_id_regex,
                "patient_id": "IPH0001",
                "sample_type": "T",
                "tumor_sample_types": tumor_sample_types,
                "expected": polars.DataFrame(
                    {
                        "vcf": ["IPH0001-D01-T01-A01_MergedSmallVariants.genome.vcf"],
                        "sample_type": ["T"],
                    }
                ),
            },
        ]

        for test_case in test_cases:
            with self.subTest(msg=test_case["name"]):
                vcf = Vcf(
                    test_case["vcf"],
                    test_case["inpred_id_regex"],
                    test_case["tumor_sample_types"],
                )
                got = vcf.row()
                assert got.equals(test_case["expected"])
