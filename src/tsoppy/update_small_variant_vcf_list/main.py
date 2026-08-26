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
from cyvcf2 import Variant

from tsoppy.general.classes import SmallVariantGenomeVcf

# Use logger that was set up in CLI
logger = logging.getLogger(__name__)


class RecurrenceVariant:
    format_fields = {
        "dragen_2.6.2.4": {"read_depth": "DP", "allele_frequency": "AF"},
        "localapp_ruo-2.2.0.12": {"read_depth": "DP", "allele_frequency": "VF"},
    }

    def __init__(self, variant: Variant):
        self.variant = variant
        self.location = f"{variant.CHROM}:{variant.POS}"
        self.ref = variant.REF
        alt_alleles = variant.ALT
        if "<NON_REF>" in alt_alleles:
            alt_alleles.remove("<NON_REF>")
        if len(alt_alleles) == 0:
            self.alt = "."
        elif len(alt_alleles) > 1:
            logger.warning(
                f"{self.location} has more than one alternative allele: {variant.ALT}. Only the first non reference will be used."
            )
        else:
            self.alt = variant.ALT[0]

    @property
    def id(self) -> str:
        return f"{self.location}:{self.ref}>{self.alt}"

    def read_depth(self, id: str) -> int:
        if self.format_fields[id]["read_depth"] in self.variant.FORMAT:
            return self.variant.format(self.format_fields[id]["read_depth"])[0][0]
        return 0


class VariantRecurrenceSummary:
    rex = r"(?P<sample_type_code>[A,N,T]):(?P<vaf_00>\d+)\+(?P<vaf_01>\d+)\+(?P<vaf_05>\d+)\+(?P<vaf_35>\d+)=(?P<vaf_any>\d+)\/(?P<call>\d+)"

    def __init__(self, summary: str):
        match = re.search(self.rex, summary)
        if not match:
            logger.error(f"{summary} is not a valid recurrence summary.")
            raise ValueError
        else:
            self.sample_type_code = match.group("sample_type_code")
            self.vaf_00 = int(match.group("vaf_00"))
            self.vaf_01 = int(match.group("vaf_01"))
            self.vaf_05 = int(match.group("vaf_05"))
            self.vaf_35 = int(match.group("vaf_35"))
            self.vaf_any = int(match.group("vaf_any"))
            self.call = int(match.group("call"))

    def __str__(self):
        return f"{self.sample_type_code}:{self.vaf_00}+{self.vaf_01}+{self.vaf_05}+{self.vaf_35}={self.vaf_any}/{self.call}"

    def update(self, vaf: float):
        if vaf < 0.01:
            self.vaf_00 += 1
        elif vaf < 0.05:
            self.vaf_01 += 1
        elif vaf < 0.35:
            self.vaf_05 += 1
        else:
            self.vaf_35 += 1
        self.vaf_any += 1

    def add_sample(self):
        self.call += 1


# my idea was to avoid the separate vcf list which we keep as of now as the same list is the header of the variant recurrence table
# this would be the main function to handle the recurrence table parsing and adding new vcfs
class VariantRecurrenceTable:
    expected_columns = [
        "variant_id",
        "tumor_recurrence_summary",
        "normal_recurrence_summary",
        "total_recurrence_summary",
    ]
    sample_vcf_rex = (
        r"^#\[sample_vcf\]\s+(?P<sample_vcf>\/\S+)\s+(?P<sample_type>[N,T])"
    )
    sample_type_dict = {"N": "normal", "T": "tumor"}
    update = False

    def __init__(self, file_path: str):
        self.file_path = Path(file_path)
        if not self.file_path.parent.exists():
            logger.error(f"{self.file_path.parent} does not exist")
            raise FileNotFoundError
        if not self.file_path.exists() or self.file_path.stat().st_size == 0:
            logger.info(
                f"{self.file_path} does not exist or is empty - will be created from scratch"
            )
            self._new_table()
        else:
            self._parse_header()
            self._parse_body()
            self.update = True

    def _new_table(self):
        """initialize empty dataframes and set info_header."""
        self.sample_vcf_header = polars.DataFrame({"sample_vcf": [], "sample_type": []})
        self.info_header = [
            '# The recurrence summary values are written in format "X:A+B+C+D=M/N"',
            '#   - X: sample type, one of "T" = tumor, "N" = normal, "A" = any',
            "#   - A: number of samples of type X, in which given variant was seen with VAF < 0.01",
            "#   - B: number of samples of type X, in which given variant was seen with 0.01 <= VAF < 0.05",
            "#   - C: number of samples of type X, in which given variant was seen with 0.05 <= VAF < 0.35",
            "#   - D: number of samples of type X, in which given variant was seen with 0.35 <= VAF",
            "#   - M: number of samples of type X, in which given variant was seen with any VAF",
            "#   - N: number of investigated samples of type X, in which given variant was callable (i.e., the variant site had coverage >= 20)",
        ]
        self.body = polars.DataFrame(dict.fromkeys(self.expected_columns, []))

    def _parse_body(self):
        """Parse table body and check if expected columns are present."""
        self.body = polars.read_csv(
            source=self.file_path, separator="\t", comment_prefix="#"
        )
        if not set(self.expected_columns).issubset(self.body.columns):
            logger.error(
                f"Expected {self.expected_columns} to be present in variant recurrence table - got {self.body.columns}"
            )
            raise ValueError
        duplicated_variants = self.body.filter(polars.col("variant_id").is_duplicated())
        duplicated_variants_set = set(duplicated_variants["variant_id"].to_list())
        if len(duplicated_variants_set) > 0:
            logger.error(
                f"Variant recurrence table contains the following duplicated variant ids: {duplicated_variants_set}."
            )
            raise ValueError
        variant_count = len(self.body)
        logger.info(f"Variant recurrence table contains {variant_count} variants.")

    def _parse_header(self):
        """Parse all header lines prefixed with '#'."""

        # Initialize dataframe dict for sample vcf lines
        dataframe_dict = {"sample_vcf": [], "sample_type": []}
        self.info_header = []
        with open(self.file_path, "r") as table:
            for line in table:
                match = re.search(self.sample_vcf_rex, line)

                # For sample vcf lines parse them into the dataframe_dict
                if match:
                    dataframe_dict["sample_vcf"].append(match.group("sample_vcf"))
                    dataframe_dict["sample_type"].append(match.group("sample_type"))

                # Other header lines are stripped of whitespace and added to the info_header
                elif line.startswith("#"):
                    self.info_header.append(line.strip())

                # We assume that header lines are always prefixed with '#' and on the very top of the file hence we break after all have been processed
                else:
                    break
        self.sample_vcf_header = polars.DataFrame(dataframe_dict)
        tumor_count = len(
            self.sample_vcf_header.filter(polars.col("sample_type") == "T")
        )
        normal_count = len(
            self.sample_vcf_header.filter(polars.col("sample_type") == "N")
        )
        any_count = len(self.sample_vcf_header)
        logger.info(
            f"Variant recurrence table contains {any_count} samples: {tumor_count} tumor and {normal_count} normal samples."
        )

    def add_vcf(self, vcf: SmallVariantGenomeVcf, min_read_depth: int = 20):
        sample_type_code = vcf.simple_sample_type_code
        sample_type = self.sample_type_dict[sample_type_code]
        if not sample_type_code in ["N", "T"]:
            logger.warning(
                f"{vcf.path} has type {sample_type_code} and cannot be added to variant recurrence table."
            )
            return

        # Initialize columns to record if sample callable at genomic location
        callable_col_name = f"{sample_type}_callable"
        self.body = self.body.with_columns(polars.lit(0).alias(callable_col_name))

        required_format_fields = self.required_format_fields[vcf.workflow_id]
        allele_frequency_field = required_format_fields["allele_frequency"]
        for variant in vcf.vcf:
            if self._variant_of_interest(
                variant, required_format_fields, min_read_depth
            ):
                var = RecurrenceVariant(variant)
                self.body = self.body.with_columns(
                    polars.when(polars.col("variant_id").str.starts_with(var.location))
                    .then(1)
                    .otherwise(polars.col(callable_col_name))
                    .alias(callable_col_name)
                )

                # Find all rows that match the variant
                rows = self.body.filter(polars.col("variant_id") == var.id)

                # Check for duplicates
                if len(rows) > 1:
                    logger.error(
                        f"Variant recurrence table contains {len(rows)} rows for {var.id}."
                    )
                    raise ValueError

                # Add row if variant not present
                elif len(rows) == 0:
                    row_dict = {
                        "variant_id": [var.id],
                        "tumor_recurrence_summary": ["T:0+0+0+0=0/0"],
                        "normal_recurrence_summary": ["N:0+0+0+0=0/0"],
                        "total_recurrence_summary": ["A:0+0+0+0=0/0"],
                        callable_col_name: [1],
                    }

                    # this should only happen if the variant is not non-ref/.
                    summary = VariantRecurrenceSummary(
                        row_dict[f"{sample_type}_recurrence_summary"][0]
                    )
                    summary.update(variant.format(allele_frequency_field)[0][0])
                    row_dict[f"{sample_type}_recurrence_summary"] = [str(summary)]
                    self.body = polars.concat([self.body, polars.DataFrame(row_dict)])
                else:
                    row = rows.row(0, named=True)
                    summary = VariantRecurrenceSummary(
                        row[f"{sample_type}_recurrence_summary"]
                    )
                    summary.update(variant.format(allele_frequency_field)[0][0])
                    self.body = self.body.with_columns(
                        polars.when(polars.col("variant_id") == variant_id)
                        .then(str(summary))
                        .otherwise(polars.col(f"{sample_type}_recurrence_summary"))
                        .alias(f"{sample_type}_recurrence_summary")
                    )

                return
        return

    def _parse_alt_allele(self, variant: Variant) -> str:
        alleles = variant.ALT
        if "<NON_REF>" in alleles:
            alleles.remove("<NON_REF>")
        if len(alleles) == 0:
            return ""
        if len(alleles) > 1:
            logger.warning(
                f"{variant.CHROM}:{variant.POS} has more than one alternate allele: {variant.ALT}. Only the first non refeference will be used."
            )
        return alleles[0]

    def _variant_of_interest(
        self,
        variant: Variant,
        required_format_fields: dict[str, str],
        min_read_depth: int,
    ) -> bool:
        """Check if variant is of interest for recurrence table."""
        if not set(required_format_fields.values()).issubset(variant.FORMAT):
            return False
        return variant.format("DP")[0][0] < min_read_depth


def Recurrence_row_valid(
    tumor: tuple[int], normal: tuple[int], any: tuple[int]
) -> bool:
    """Check if sample numbers for the recurrence of a variant is consistent across one row."""
    if len(tumor) != len(normal) or len(tumor) != len(any):
        return False
    for i in range(len(tumor)):
        if tumor[i] + normal[i] != any[i]:
            return False
    return True


# from here on it is old code and should be disregarded
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
        return self.vcfs == other.vcfs

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
        return self.vcf == other.vcf

    def row(self):
        """
        Return small variant VCF list row.
        """
        return polars.DataFrame({"vcf": [self.vcf], "sample_type": [self.sample_type]})
