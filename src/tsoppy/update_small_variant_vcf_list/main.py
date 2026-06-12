"""
This module defines the classes 'VcfList' and 'Vcf'.
'VcfList' takes a directory holding TSO500 results, a glob to identify small variant vcf files, the currect small variant vcf list,
a regular expression matching InPreD IDs, a set of tumor sample types and the path to the new small variant vcf list.
'VcfList' has a method to update the current small variant vcf list with vcfs found in the TSO500 results directory.
'Vcf' is defined by a path to a small variant vcf and a set of tumor sample types.
'Vcf' provides a method to create a new row for a polars dataframe.
"""

import glob
import logging
import re
from datetime import datetime
from pathlib import Path

import polars

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
        vcfs (dict): Small variant VCF(s) located in TSO500 results directory.
    """

    vcf_list_columns = {"vcf": polars.String, "sample_type": polars.String}

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
        self.dataframe = polars.DataFrame(schema=self.vcf_list_columns)

        # Try reading small variant VCF list or start from scratch
        if vcf_list:
            try:
                self.dataframe = polars.read_csv(
                    source=vcf_list,
                    separator="\t",
                    schema=self.vcf_list_columns,
                    ignore_errors=True,
                    has_header=False,
                    raise_if_empty=False,
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
        if not self.dataframe.equals(other.dataframe):
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
            if small_variant_vcf.vcf in self.dataframe["vcf"].to_list():
                logger.warning(f"{vcf} is already in small variant VCF list, skipping.")
                continue

            # Exclude control samples
            if small_variant_vcf.patient_id.startswith("IPC"):
                logger.warning(
                    f"{small_variant_vcf.patient_id} is a control sample, skipping."
                )
                continue

            # Add vcf to list
            self.dataframe = polars.concat([self.dataframe, small_variant_vcf.row()])

            # Check if new patient ID is represented multiple times
            patient_sample_count = (
                self.dataframe["vcf"] == small_variant_vcf.patient_id
            ).sum()
            if patient_sample_count > 1:
                logger.warning(
                    f"patient {small_variant_vcf.patient_id} has {patient_sample_count} vcf(s) in the small variant VCF list."
                )

        # Write updated small variant VCF list to file
        self.dataframe.unique().write_csv(
            file=self.output, separator="\t", include_header=False
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
        return polars.DataFrame({"vcf": [self.vcf], "sample_type": [self.sample_type]})
