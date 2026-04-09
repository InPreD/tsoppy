"""
This module contains the code for the `update_small_variant_vcf_list` command.
The command takes two arguments, `results_dir`, which is a string that specifies the directory where the results of the latest TSO500 run are stored.
"""

import glob
import logging
import re
from datetime import datetime
from pathlib import Path

import pandas

# Use logger that was set up in CLI
logger = logging.getLogger(__name__)


class InvalidSampleType(Exception):
    """
    Exception if sample type is not valid.
    """

    def __init__(self, msg="sample type is not valid"):
        self.msg = msg
        super().__init__(self.msg)

    def __str__(self):
        return self.msg


class VcfList:
    """
    Represents small variant VCF list.

    Attributes:
        dataframe (Dataframe): Dataframe representing the current version of small variant VCF list.
        inpred_id_regex (str): Regular expression matching InPreD IDs.
        output (str): Path to updated version of small variant VCF list.
        tumor_sample_types (set[str]): Single letter codes representing a tumor sample.
        vcf_list_columns (list[str]): List of dataframe column names.
        vcfs (list[str]): Small variant VCF(s) located in TSO500 results directory.
    """

    vcf_list_columns = ["vcf", "sample_type"]

    def __init__(
        self,
        results_dir: Path,
        glob_pattern: str,
        vcf_list: Path | None,
        inpred_id_regex: str,
        tumor_sample_types: str,
        output: str,
    ):
        """
        Create new instance of SmallVariantVcfList.
        """
        self.vcfs = glob.glob(f"{results_dir}/{glob_pattern}")
        self.inpred_id_regex = rf"{inpred_id_regex}"
        self.tumor_sample_types = set(tumor_sample_types.split(","))
        self.dataframe = pandas.DataFrame(columns=self.vcf_list_columns)

        # Try reading small variant VCF list or start from scratch
        if vcf_list:
            try:
                self.dataframe = pandas.read_csv(
                    vcf_list, sep="\t", names=self.vcf_list_columns, on_bad_lines="warn"
                )
            except FileNotFoundError:
                logger.warning(
                    f"{vcf_list} not found, creating new small variant VCF list."
                )
        else:
            logger.info("no small variant VCF list specified, creating new one.")

        # Replace placeholder with actual date
        if "<YYYYMMDD>" in output:
            now = datetime.now()
            self.output = output.replace("<YYYYMMDD>", now.strftime("%Y%m%d"))
        else:
            self.output = output

    def __eq__(self, other):
        """
        Compare to other class instance.
        """
        if not isinstance(other, VcfList):
            return NotImplemented
        if self.dataframe != other.dataframe:
            return False
        if self.inpred_id_regex != other.inpred_id_regex:
            return False
        if self.output != other.output:
            return False
        if self.tumor_sample_types != other.tumor_sample_types:
            return False
        return self.vcfs != other.vcfs

    def update(self):
        """
        Add VCF(s) from results directory to small variant VCF list.
        """

        # Loop over all small variant VCFs
        for vcf in self.vcfs:
            # Try to create vcf class instance
            try:
                small_variant_vcf = Vcf(
                    vcf, self.inpred_id_regex, self.tumor_sample_types
                )
            except AttributeError:
                logger.warning(
                    f"could not parse InPreD ID from {small_variant_vcf.vcf}, skipping."
                )
                continue
            except InvalidSampleType:
                logger.warning(
                    f"{small_variant_vcf.vcf} has sample type {small_variant_vcf.sample_type} which is not {self.tumor_sample_types} or N(ormal), skipping."
                )
                continue

            # Avoid duplication
            if small_variant_vcf.vcf in self.dataframe["vcf"].values:
                logger.warning(f"{vcf} is already in small variant VCF list, skipping.")
                continue

            # Exclude control samples
            if small_variant_vcf.patient_id.startswith("IPC"):
                logger.warning(
                    f"{small_variant_vcf.patient_id} is a control sample, skipping."
                )
                continue

            # Add vcf to list
            self.dataframe.loc[len(self.dataframe)] = small_variant_vcf.row()

            # Check if new patient ID is represented multiple times
            patient_sample_count = (
                self.dataframe["vcf"].str.contains(small_variant_vcf.patient_id).sum()
            )
            if patient_sample_count > 1:
                logger.warning(
                    f"patient {small_variant_vcf.patient_id} has {patient_sample_count} vcf(s) in the small variant VCF list."
                )

        # Write updated small variant VCF list to file
        self.dataframe.drop_duplicates().to_csv(
            self.output, sep="\t", header=False, index=False
        )


class Vcf:
    """
    Represents small variant VCF.

    Attributes:
        patient_id (str): ID of patient that the VCF belongs to.
        sample_type (str): Single letter code representing type of sample, e.g. T = tumor.
        vcf (str): Path to VCF file.
    """

    def __init__(self, vcf: str, inpred_id_regex: str, tumor_sample_types: set):
        """
        Create new instance of SmallVariantVcf.
        """
        self.vcf = vcf

        # Parse InPreD ID to get patient ID and sample type
        match = re.search(inpred_id_regex, self.vcf)
        try:
            self.patient_id = match.group("patient_id")
            self.sample_type = match.group("sample_type")
        except AttributeError:
            raise AttributeError

        # Validate sample type is N(ormal) or included in tumor_sample_types
        if self.sample_type != "N":
            if self.sample_type not in tumor_sample_types:
                raise InvalidSampleType
            else:
                # Reset any sample type in tumor_sample_types with T
                logger.warning(
                    f"sample type code {self.sample_type} for {self.vcf} will be replaced with T"
                )
                self.sample_type = "T"

    def __eq__(self, other):
        """
        Compare to other class instance.
        """
        if not isinstance(other, Vcf):
            return NotImplemented
        if self.patient_id != other.patient_id:
            return False
        if self.sample_type != other.sample_type:
            return False
        return self.vcf != other.vcf

    def row(self):
        """
        Return small variant VCF list row.
        """
        return [self.vcf, self.sample_type]
