"""
update small variant vcf list subpackage main module unit tests.
"""
from os import path
import os
import filecmp
import unittest
from tsoppy.update_small_variant_vcf_list.main import VcfList, Vcf

# Define path to test data - cannot be absolute due to different paths locally and in CI
test_data_dir = "tests/test_data/update_small_variant_vcf_list_main"

# test constants
glob_pattern = "**/Results/**/*_MergedSmallVariants.genome.vcf"
tumor_sample_types = "C,D,d,L,M,P,p,R,r,T,X"
inpred_id_regex = "(?P<patient_id>\\D{3}\\d{4})-\\D\\d{2}-(?P<sample_type>\\D)\\d{2}-\\D\\d{2}.*.vcf$"


class TestVcfList(unittest.TestCase):
    def test_update(self):
        test_cases = [
            {
                "name": "successfully update small variant vcf list",
                "results_dir": path.join(test_data_dir, "successfully_update_small_variant_vcf_list"),
                "glob_pattern": glob_pattern,
                "vcf_list": path.join(test_data_dir, "successfully_update_small_variant_vcf_list/TSO500_vcf_list.tsv"),
                "inpred_id_regex": inpred_id_regex,
                "tumor_sample_types": tumor_sample_types,
                "output": path.join(test_data_dir, "successfully_update_small_variant_vcf_list/TSO500_vcf_list_updated.tsv"),
                "expected": path.join(test_data_dir, "successfully_update_small_variant_vcf_list/TSO500_vcf_list_expected.tsv"),
            },
            {
                "name": "create new small variant vcf list",
                "results_dir": path.join(test_data_dir, "create_new_small_variant_vcf_list"),
                "glob_pattern": glob_pattern,
                "vcf_list": None,
                "inpred_id_regex": inpred_id_regex,
                "tumor_sample_types": tumor_sample_types,
                "output": path.join(test_data_dir, "create_new_small_variant_vcf_list/TSO500_vcf_list_updated.tsv"),
                "expected": path.join(test_data_dir, "create_new_small_variant_vcf_list/TSO500_vcf_list_expected.tsv"),
            },
            {
                "name": "small variant vcf list does not exist",
                "results_dir": path.join(test_data_dir, "small_variant_vcf_list_does_not_exist"),
                "glob_pattern": glob_pattern,
                "vcf_list": path.join(test_data_dir, "small_variant_vcf_list_does_not_exist/TSO500_vcf_list.tsv"),
                "inpred_id_regex": inpred_id_regex,
                "tumor_sample_types": tumor_sample_types,
                "output": path.join(test_data_dir, "small_variant_vcf_list_does_not_exist/TSO500_vcf_list_updated.tsv"),
                "expected": path.join(test_data_dir, "small_variant_vcf_list_does_not_exist/TSO500_vcf_list_expected.tsv"),
            },
            {
                "name": "skip existing vcf",
                "results_dir": path.join(test_data_dir, "skip_existing_vcf"),
                "glob_pattern": glob_pattern,
                "vcf_list": path.join(test_data_dir, "skip_existing_vcf/TSO500_vcf_list.tsv"),
                "inpred_id_regex": inpred_id_regex,
                "tumor_sample_types": tumor_sample_types,
                "output": path.join(test_data_dir, "skip_existing_vcf/TSO500_vcf_list_updated.tsv"),
                "expected": path.join(test_data_dir, "skip_existing_vcf/TSO500_vcf_list_expected.tsv"),
            },
            {
                "name": "inpred id not parsable",
                "results_dir": path.join(test_data_dir, "inpred_id_not_parsable"),
                "glob_pattern": glob_pattern,
                "vcf_list": path.join(test_data_dir, "inpred_id_not_parsable/TSO500_vcf_list.tsv"),
                "inpred_id_regex": inpred_id_regex,
                "tumor_sample_types": tumor_sample_types,
                "output": path.join(test_data_dir, "inpred_id_not_parsable/TSO500_vcf_list_updated.tsv"),
                "expected": path.join(test_data_dir, "inpred_id_not_parsable/TSO500_vcf_list_expected.tsv"),
            },
        ]

        for test_case in test_cases:
            with self.subTest(msg=test_case["name"]):
                got = VcfList(
                    test_case["results_dir"], test_case["glob_pattern"], test_case["vcf_list"], test_case["inpred_id_regex"], test_case["tumor_sample_types"], test_case["output"])
                got.update()
                assert filecmp.cmp(test_case["output"], test_case["expected"])
                os.remove(test_case["output"])


class TestVcf(unittest.TestCase):
    def test_init(self):
        test_cases = [
            {
                "name": "include sample",
                "vcf": "IPH0001-01-T01-01_MergedSmallVariants.genome.vcf",
                "patient_id": "IPH0001",
                "sample_type": "T",
                "tumor_sample_types": tumor_sample_types,
                "expected": True,
            },
            {
                "name": "sample is control",
                "vcf": "IPC0001-01-T01-01_MergedSmallVariants.genome.vcf",
                "patient_id": "IPC0001",
                "sample_type": "T",
                "tumor_sample_types": tumor_sample_types,
                "expected": False,
            },
            {
                "name": "sample is neither tumor nor normal",
                "vcf": "IPH0001-01-A01-01_MergedSmallVariants.genome.vcf",
                "patient_id": "IPH0001",
                "sample_type": "A",
                "tumor_sample_types": tumor_sample_types,
                "expected": False,
            },
        ]

        for test_case in test_cases:
            with self.subTest(msg=test_case["name"]):
                got = Vcf(
                    test_case["vcf"], test_case["patient_id"], test_case["sample_type"], test_case["tumor_sample_types"])
                assert got.include == test_case["expected"]

    def test_row(self):
        test_cases = [
            {
                "name": "successfully return row",
                "vcf": "IPH0001-01-T01-01_MergedSmallVariants.genome.vcf",
                "patient_id": "IPH0001",
                "sample_type": "T",
                "tumor_sample_types": tumor_sample_types,
                "expected": ["IPH0001-01-T01-01_MergedSmallVariants.genome.vcf", "T"],
            },
        ]

        for test_case in test_cases:
            with self.subTest(msg=test_case["name"]):
                vcf = Vcf(
                    test_case["vcf"], test_case["patient_id"], test_case["sample_type"], test_case["tumor_sample_types"])
                got = vcf.row()
                assert got == test_case["expected"]
