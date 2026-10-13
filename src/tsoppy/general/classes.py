import gzip
import logging
import os
import re
from pathlib import Path

import cyvcf2
import msgspec
import polars

from tsoppy.general.file_parser import Parse_section_csv, Parse_section_tsv

# Use logger that was set up in CLI
logger = logging.getLogger(__name__)


class WorkflowConfig(msgspec.Struct):
    """Config class for workflow output file path format strings.

    Attributes:
        metrics_output_tsv: Paths to metrics output tsv (dict[str, str])
        samplesheet: Paths to samplesheet (dict[str, str])
        small_variant_genome_vcf: Paths to small variant genome vcf (dict[str, str])
        tmb_trace_tsv: Paths to tmb trace tsv (dict[str, str])
        variants_annotated_json: Paths to variants annotated json (dict[str, str])
    """

    metrics_output_tsv: dict[str, str]
    samplesheet: dict[str, str]
    small_variant_genome_vcf: dict[str, str]
    tmb_trace_tsv: dict[str, str]
    variants_annotated_json: dict[str, str]

    def __eq__(self, other):
        """Assess if two instances of this class are equal."""
        if not isinstance(other, WorkflowConfig):
            return False
        attr_to_compare = [
            "metrics_output_tsv",
            "samplesheet",
            "small_variant_genome_vcf",
            "tmb_trace_tsv",
            "variants_annotated_json",
        ]
        return all(vars(self).get(k) == vars(other).get(k) for k in attr_to_compare)


class InPredNomenclatureCode(msgspec.Struct):
    """Class for inpred nomenclature code.

    Attributes:
        code: Letters representing the nomenclature code (str)
        choices: Possible choices for nomenclature code (dict[str, str])
    """

    code: str
    choices: dict[str, str]

    def __eq__(self, other):
        """Assess if two instances of this class are equal."""
        if not isinstance(other, InPredNomenclatureCode):
            return False
        attr_to_compare = ["code", "choices"]
        return all(vars(self).get(k) == vars(other).get(k) for k in attr_to_compare)


class InPredNomenclature(msgspec.Struct):
    """Class for inpred nomenclature.

    Attributes:
        format: sample id format (str)
        project_code: Code representing project (InPredNomenclatureCode)
        patient_code: Code representing patient (dict[str, str])
        nucleic_acid_input_type_code: Code representing nucleic acid input type (InPredNomenclatureCode)
        assay_type_code: Code representing assay type (InPredNomenclatureCode)
        sample_type_code: Code representing sample type (InPredNomenclatureCode)
        library_preparation_attempt: Code representing library preparation (InPredNomenclatureCode)
        biological_replicate: Code representing biological replicate (dict[str, str])
        sample_material: Code representing sample material (InPredNomenclatureCode)
        tumor_site: Code representing tumor site (InPredNomenclatureCode)
    """

    format: str
    project_code: InPredNomenclatureCode
    patient_code: dict[str, str]
    nucleic_acid_input_type_code: InPredNomenclatureCode
    assay_type_code: InPredNomenclatureCode
    sample_type_code: InPredNomenclatureCode
    library_preparation_attempt: InPredNomenclatureCode
    biological_replicate: dict[str, str]
    sample_material: InPredNomenclatureCode
    tumor_site: InPredNomenclatureCode

    def __eq__(self, other):
        """Assess if two instances of this class are equal."""
        if not isinstance(other, InPredNomenclature):
            return False
        attr_to_compare = [
            "format",
            "project_code",
            "patient_code",
            "nucleic_acid_input_type_code",
            "assay_type_code",
            "sample_type_code",
            "library_preparation_attempt",
            "biological_replicate",
            "sample_material",
            "tumor_site",
        ]
        return all(vars(self).get(k) == vars(other).get(k) for k in attr_to_compare)


class WorkflowOutput:
    """Base class for outputs produced by different workflows (e.g. dragen/localapp).

    Attributes:
        config: Configuration (WorkflowConfig)
        nomenclature: InPred nomenclature (InPredNomenclature)
        root: Root path (Path)
        samples: Data section of samplesheet (polars.DataFrame)
        workflow_type: Detected workflow type (str)
        workflow_version: Detected workflow version (str)
    """

    sample_id_rex = r"^(?P<project>\w{3})(?P<patient>\d{4})-(?P<input>\w)(?P<assay>\d{2})-(?P<sample>\w)(?P<preparation>\d)(?P<replicate>\d)-(?P<material>\w)(?P<tumor>\d{2})"

    def __init__(
        self,
        config_yaml: str | Path,
        inpred_nomenclature: str | Path,
        root_path: str | Path,
    ):
        """Initialize WorkflowOutput."""
        self.root = Path(root_path)
        with open(config_yaml, "r") as yaml_file:
            self.config = msgspec.yaml.decode(yaml_file.read(), type=WorkflowConfig)
        with open(inpred_nomenclature, "r") as yaml_file:
            self.nomenclature = msgspec.yaml.decode(
                yaml_file.read(), type=InPredNomenclature
            )

        self._detect_type_and_version()
        self._parse_samples()

    def __eq__(self, other):
        """Assess if two instances of this class are equal."""
        if not isinstance(other, WorkflowOutput):
            return False
        if self.config != other.config:
            return False
        if self.samples.equals(other.samples):
            return False
        attr_to_compare = ["root", "workflow_type", "workflow_version"]
        return all(vars(self).get(k) == vars(other).get(k) for k in attr_to_compare)

    def _detect_type_and_version(self):
        """Detect which workflow type and version is present based on information in MetricsOutput.tsv."""

        # Get all values for MetricsOutput.tsv paths and check if they are the same
        info_src = list(self.config.metrics_output_tsv.values())
        if len(set(info_src)) != 1:
            logger.error(
                f"Got {info_src} but need exactly one file to detect workflow id."
            )
            raise ValueError

        # Parse MetricsOutput.tsv
        headers, sections = Parse_section_tsv(
            os.path.join(self.root, info_src[0]), ["Header"]
        )

        # Check if DRAGEN is part of the header and assume the data is localapp if not
        if "DRAGEN" in headers[0]:
            self.workflow_type = "dragen"
        else:
            self.workflow_type = "localapp"

        # Set workflow version from Header section
        self.workflow_version = sections["Header"].item(
            row=0, column="Workflow Version"
        )

    def _parse_samples(self):
        """Parse samples from samplesheet."""
        samplesheet_path = os.path.join(
            self.root, self.config.samplesheet[self.workflow_id]
        )

        # Parse SampleSheet.csv
        _, sections = Parse_section_csv(samplesheet_path, [])

        # Find and check that there is only one Data section present
        data_sections = [section for section in list(sections) if "Data" in section]
        if len(data_sections) != 1:
            logger.error(
                f"Expected one data section in samplesheet {samplesheet_path} - got {len(data_sections)}."
            )
            raise KeyError

        # Ensure that the table contains a Sample_ID column
        if "Sample_ID" not in sections[data_sections[0]]:
            logger.error(
                f"Sample_ID column is missing from samplesheet {samplesheet_path}."
            )
            raise KeyError

        # Store the complete table
        self.samples = sections[data_sections[0]]

    def _translate_code(self, code: str, choice: str) -> str:
        code_dict = getattr(self.nomenclature, code)
        if not hasattr(code_dict, "choices"):
            return choice
        else:
            choices = list(code_dict.choices.keys())
            if choice not in choices:
                logger.warning(f"{choice} is not valid - expected any: {choices}")
                return choice
            else:
                return code_dict.choices[choice]

    def sample_exists(self, sample_id: str) -> bool:
        """Check if sample exists."""
        return sample_id in self.sample_list

    @property
    def sample_list(self) -> list[str]:
        """Return all samples present in samplsheet."""
        return self.samples["Sample_ID"].to_list()

    def sample_meta(self, sample_id: str) -> dict[str, str] | None:
        """Parse meta information from sample id and return as dict"""
        match = re.search(self.sample_id_rex, sample_id)
        if not match:
            logger.warning(
                "Sample {sample_id} does not follow the InPreD sample ID nomenclature."
            )
            return None
        return {
            "project": self._translate_code("project_code", match.group("project")),
            "patient": self._translate_code("patient_code", match.group("patient")),
            "nucleic_acid_input_type": self._translate_code(
                "nucleic_acid_input_type_code", match.group("input")
            ),
            "assay_type": self._translate_code("assay_type_code", match.group("assay")),
            "sample_type": self._translate_code(
                "sample_type_code", match.group("sample")
            ),
            "library_preparation_attempt": self._translate_code(
                "library_preparation_attempt", match.group("preparation")
            ),
            "biological_replicate": self._translate_code(
                "biological_replicate", match.group("replicate")
            ),
            "sample_material": self._translate_code(
                "sample_material", match.group("material")
            ),
            "tumor_site": self._translate_code("tumor_site", match.group("tumor")),
        }

    @property
    def workflow_id(self) -> str:
        """Return combined string for workflow type and version."""
        return f"{self.workflow_type}_{self.workflow_version}"


class SmallVariantGenomeVcf(WorkflowOutput):
    """Input class for small variant genome VCF files produced by different workflows.

    Attributes:
        header_rex: Regex to parse top level key value pairs in VCF object header (str)
        header_subrex: Regex to parse sub level key value pairs in VCF object header (str)
        path: Path to vcf (Path)
        sample_id: Sample identifier (str)
        vcf: Parsed VCF object (cyvcf2.VCF)
        merged_table: Merged small variant data (polars.DataFrame)
    """

    # Define regex patterns for parsing VCF header lines and storing in a dictionary
    header_dict = {}
    # Top level Regex (i.e. '##key=value')
    header_rex = r"##(?P<key>\w+)=(?P<value>.+)"
    # Sub level Regex (parsing the top-level "value" field)
    header_subrex = r"(\w+)=([\w\d]+|\".+\")?"

    def __init__(
        self,
        config_yaml: str | Path,
        inpred_nomenclature: str | Path,
        root_path: str | Path,
        sample_id: str,
    ):
        """Initialize SmallVariantGenomeVcf"""
        super().__init__(config_yaml, inpred_nomenclature, root_path)
        self.sample_id = sample_id
        self._parse()

    @classmethod
    def create(cls, workflow_output: WorkflowOutput, sample_id: str):
        """Create SmallVariantGenomeVcf from existing WorkflowOutput"""
        if not workflow_output.sample_exists(sample_id):
            logger.error(f"Sample {sample_id} does not exist")
            raise ValueError
        else:
            obj = cls.__new__(cls)
            obj.__dict__.update(workflow_output.__dict__)
            obj.sample_id = sample_id
            obj._parse()
            return obj

    def merge(
        self,
        workflow_output: WorkflowOutput,
        output_path: str | Path | None = None,
    ) -> None | polars.DataFrame:
        """Merge VCF enriched with TMB and Nirvana annotations and either save to file or return polars dataframe."""

        tmb_obj = TmbTraceTsv.create(workflow_output, self.sample_id)
        json_obj = VariantsAnnotatedJson.create(workflow_output, self.sample_id)

        # Write to output path if specified, otherwise return merged table as polars dataframe
        if output_path:
            self._write_vcf(tmb_obj, json_obj, output_path)
        else:
            self.merged_table = self._to_dataframe(tmb_obj, json_obj)
            return self.merged_table

    def _write_vcf(
        self,
        tmb_obj: TmbTraceTsv,
        json_obj: VariantsAnnotatedJson,
        output_path: str | Path,
    ):
        """Write a merged VCF with TMB and Nirvana annotations."""

        # Specify new INFO fields to be added to the VCF header
        info_headers = (
            {
                "ID": "variant_ID",
                "Number": "1",
                "Type": "String",
                "Description": "Variant identifier, composed of CHROM, POS, REF and ALT fields",
            },
            {
                "ID": "variant_type",
                "Number": "1",
                "Type": "String",
                "Description": "Type of genomic change: 'MNV', 'SNV', 'insertion', 'deletion' or 'indel'",
            },
            {
                "ID": "overlapping_genes",
                "Number": ".",
                "Type": "String",
                "Description": "List of genes overlapping the variant site, as reported by Nirvana in form of HGNC gene symbols",
            },
            {
                "ID": "Illumina_variant_class",
                "Number": "1",
                "Type": "String",
                "Description": "Variant classification provided by initial Illumina pipeline analysis: 'Germline_DB', 'Germline_Proxi', 'Somatic', 'Blacklist', 'VCF_filtered'",
            },
        )

        for info in info_headers:
            self.vcf.add_info_to_header(info)

        # Create dictionaries for TMB and Nirvana annotations for mapping overlaps
        tmb_by_variant = tmb_obj.variant_dict
        tmb_df = tmb_obj.class_table
        json_by_variant = json_obj.variant_dict

        # Ensure output path is a Path object for cyvcf2.Writer
        out_path = Path(output_path)
        writer = cyvcf2.Writer(out_path, self.vcf)
        for variant in self.vcf:
            # Skip variants with no alternate allele or only <NON_REF> alternate allele
            if not variant.ALT or variant.ALT == ["<NON_REF>"]:
                continue

            # Generate variant ID and add annotations to INFO fields
            variant_id = get_variant_id(
                variant.CHROM, variant.POS, variant.REF, variant.ALT[0]
            )
            variant.INFO["variant_ID"] = variant_id

            # Get variant type and overlapping genes from Nirvana annotations
            variant_type, genes = (
                json_by_variant[variant_id]
                if variant_id in json_by_variant
                else ("", [])
            )
            variant.INFO["variant_type"] = variant_type
            variant.INFO["overlapping_genes"] = ",".join(genes)

            # Determine Illumina variant class based on TMB annotations and VCF filter status
            if variant_id in tmb_by_variant:
                variant_class = (
                    tmb_df.filter(polars.col("variant_ID") == variant_id)
                    .select("Class")
                    .item()
                )
            elif variant.FILTER == "excluded_regions" or variant.FILTER == "Blacklist":
                variant_class = "Blacklist"
            else:
                variant_class = "VCF_filtered"

            variant.INFO["Illumina_variant_class"] = variant_class

            writer.write_record(variant)

        writer.close()
        self.vcf.close()
        logger.info(f"Merged VCF file saved to '{out_path}'")

    def _to_dataframe(
        self, tmb_obj: TmbTraceTsv, json_obj: VariantsAnnotatedJson
    ) -> polars.DataFrame:
        """Return merged variant data as a Polars dataframe."""

        # Create dictionary for storing variant information from VCF, TMB, and Nirvana annotations, to be converted to a Polars dataframe
        json_by_variant = json_obj.variant_dict
        var_stats = ["AD", "DP", "AF" if self.workflow_type == "dragen" else "VF"]
        cols = (
            ["variant_ID", "variant_type", "overlapping_genes"] + var_stats + ["filter"]
        )
        var_dict = {c: [] for c in cols}

        # Iterate over variants in the VCF file
        for variant in self.vcf:
            # Skip variants with no alternate allele or only <NON_REF> alternate allele
            if not variant.ALT or variant.ALT == ["<NON_REF>"]:
                continue

            # Generate variant ID and add annotations to INFO fields
            variant_id = get_variant_id(
                variant.CHROM, variant.POS, variant.REF, variant.ALT[0]
            )

            # Get variant type and overlapping genes from Nirvana annotations
            variant_type, genes = (
                json_by_variant[variant_id]
                if variant_id in json_by_variant
                else ("", [])
            )

            # Add variant information to the dictionary for conversion to Polars dataframe
            var_dict["variant_ID"].append(variant_id)
            var_dict["variant_type"].append(variant_type)
            var_dict["overlapping_genes"].append(genes)
            var_dict["filter"].append(variant.FILTER)
            for stat in var_stats:
                var_dict[stat].append(variant.format(stat)[0])

        # Convert raw variant information dictionary to a Polars dataframe
        df = polars.DataFrame(data=var_dict)

        # Check for single-element list columns (e.g. DP-column), remove list structure if found
        explode_cols = [
            name
            for name, dtype in df.schema.items()
            if isinstance(dtype, polars.Array) and dtype.shape == (1,)
        ]
        if explode_cols:
            df = df.explode(explode_cols, empty_as_null=True)

        # Join the VCF dataframe with the TMB class table on variant_ID
        joined_df = df.join(
            tmb_obj.class_table,
            left_on="variant_ID",
            right_on="variant_ID",
            how="left",
            maintain_order="left",
        )

        # Define the column order for the final merged dataframe
        col_order = [
            "variant_ID",
            "variant_type",
            "Illumina_variant_class",
            "overlapping_genes",
        ] + var_stats

        if self.workflow_type == "dragen":
            full_df = (
                # Fill null entries in 'filter' column with TMB class (all VCF PASS-variants)
                joined_df.with_columns(
                    polars.col("filter")
                    .fill_null(polars.col("Class"))
                    # Replace 'excluded_regions' with 'Blacklist'
                    .replace({"excluded_regions": "Blacklist"}),
                )
                # Replace any remaining filters with 'VCF_filtered'
                .with_columns(
                    polars.when(polars.col("filter").str.contains(r"^[a-z]"))
                    .then(polars.lit("VCF_filtered"))
                    .otherwise(polars.col("filter"))
                    .alias("Illumina_variant_class")
                )
                .select(col_order)
            )
        else:
            full_df = (
                # Save 'Blacklist' and 'VCF_filtered' in the 'filter' column
                joined_df.with_columns(
                    polars.when(polars.col("filter") == "Blacklist")
                    .then(polars.lit("Blacklist"))
                    .when(polars.col("filter").is_not_null())
                    .then(polars.lit("VCF_filtered"))
                    .otherwise(polars.col("filter"))
                    .alias("filter")
                )
                # Replace null-values in 'filter' column with TMB class (all VCF PASS-variants) and rename column
                .with_columns(polars.col("filter").fill_null(polars.col("Class")))
                .rename({"filter": "Illumina_variant_class"})
                .select(col_order)
                # Rename VF to AF for consistency with dragen output
                .rename({"VF": "AF"})
            )

        return full_df

    def _parse(self):
        """Parse the small variant genome VCF file"""
        fmt = self.config.small_variant_genome_vcf[self.workflow_id]
        self.path = Path(
            os.path.join(self.root, fmt.format(self.sample_id, self.sample_id))
        )
        if not self.path.is_file():
            logging.error(
                f"Small variant genome VCF missing: File {self.path} does not exist."
            )
            raise FileNotFoundError
        self.vcf = cyvcf2.VCF(self.path)
        self._parse_header()

    def _parse_header(self):
        """Parse VCF header into dict"""
        for match in re.finditer(self.header_rex, self.vcf.raw_header):
            # Use named capture groups to split line into key and value (i.e. '##key=value')
            key = match.groupdict()["key"]
            value = match.groupdict()["value"]

            # Handle command line arguments separately, as they can occur multiple times (specific for LocalApp as for now).
            if key.lower().endswith("_cmdline"):
                # Store multiple command line arguments as a list
                self.header_dict.setdefault(key, []).append(value)
                continue

            # Check if value containes additional key value pairs in the 'value' field (i.e. '##key=<subkey1=subvalue1,subkey2=subvalue2,...>')
            if not bool(re.search(self.header_subrex, value)):
                self.header_dict[key] = value
            else:
                # Parse each 'key=value' pairs into dict (e.g {'ID': 'DP', 'Number': '1', 'Type': 'Integer', 'Description': 'Approximate read depth'})
                item = dict()
                for match2 in re.finditer(self.header_subrex, value):
                    g = match2.groups()
                    item[g[0]] = g[1]

                # Use ID value as subkey if it exists, else store the key value pairs in a list under the top level key
                if "ID" not in item:
                    if key not in self.header_dict:
                        self.header_dict[key] = list()
                    self.header_dict[key].append(item)
                else:
                    if key not in self.header_dict:
                        self.header_dict[key] = dict()
                    subkey = item.pop("ID")
                    self.header_dict[key][subkey] = item


class TmbTraceTsv(WorkflowOutput):
    """Input class for TMB trace files produced by different workflows.

    Attributes:
        path: Path to vcf (Path)
        table: Parsed rows of the TMB trace file (polars.DataFrame)
        sample_id: Sample identifier (str)
        variant_dict: Variant ID mapped to parsed table (dict)
        class_table: Two-column table with Variant ID and Illumina variant class
    """

    def __init__(
        self,
        config_yaml: str | Path,
        inpred_nomenclature: str | Path,
        root_path: str | Path,
        sample_id: str,
    ):
        """Initialize TmTraceTsv."""
        super().__init__(config_yaml, inpred_nomenclature, root_path)
        self.sample_id = sample_id
        self._parse()

    @classmethod
    def create(cls, workflow_output: WorkflowOutput, sample_id: str):
        """Create TmbTraceTsv from existing WorkflowOutput."""
        if not workflow_output.sample_exists(sample_id):
            logger.error(f"Sample {sample_id} does not exist")
            raise ValueError
        else:
            obj = cls.__new__(cls)
            obj.__dict__.update(workflow_output.__dict__)
            obj.sample_id = sample_id
            obj._parse()
            return obj

    def _parse(self):
        """Parse the TMB trace tsv."""
        fmt = self.config.tmb_trace_tsv[self.workflow_id]
        self.path = Path(
            os.path.join(self.root, fmt.format(self.sample_id, self.sample_id))
        )
        if not self.path.is_file():
            logging.error(
                f"Small variant genome VCF missing: File {self.path} does not exist."
            )
            raise FileNotFoundError
        self.table = polars.read_csv(self.path, separator="\t")
        self.class_table = self._get_variant_class_table()
        self.variant_dict = {
            get_variant_id(
                row["Chromosome"], row["Position"], row["RefCall"], row["AltCall"]
            ): row
            for row in self.table.iter_rows(named=True)
        }

    def _get_variant_class_table(self):
        """Create a two-column table with Variant ID and Illumina variant class from the TMB trace file."""
        # Create a new column 'variant_ID' by concatenating Chromosome, Position, RefCall, and AltCall
        tmb_df = self.table.with_columns(
            (
                polars.concat_str(["Chromosome", "Position", "RefCall"], separator=":")
                + polars.lit(">")
                + polars.col("AltCall").cast(polars.String)
            ).alias("variant_ID")
        )
        # Determine Illumina variant class based on workflow type and TMB annotations
        if self.workflow_type == "localapp":
            tmb_df = tmb_df.with_columns(
                polars.when(polars.col("GermlineFilterDatabase"))
                .then(polars.lit("Germline_DB"))
                .when(polars.col("GermlineFilterProxi"))
                .then(polars.lit("Germline_Proxi"))
                .otherwise(polars.lit("Somatic"))
                .alias("Status")
            )
        return tmb_df.select(["variant_ID", "Status"]).rename({"Status": "Class"})


class VariantsAnnotatedJson(WorkflowOutput):
    """Input class for annotated JSON files produced by different workflows.

    Attributes:
        path: Path to vcf (Path)
        data: Parsed JSON data (dict)
        sample_id: Sample identifier (str)
        variant_dict: Variant ID mapped to type and overlapping genes (dict)
    """

    def __init__(
        self,
        config_yaml: str | Path,
        inpred_nomenclature: str | Path,
        root_path: str | Path,
        sample_id: str,
    ):
        """Initialize VariantsAnnotatedJson."""
        super().__init__(config_yaml, inpred_nomenclature, root_path)
        self.sample_id = sample_id
        self._parse()

    @classmethod
    def create(cls, workflow_output: WorkflowOutput, sample_id: str):
        """Create VariantsAnnotatedJson from existing WorkflowOutput."""
        if not workflow_output.sample_exists(sample_id):
            logger.error(f"Sample {sample_id} does not exist")
            raise ValueError
        else:
            obj = cls.__new__(cls)
            obj.__dict__.update(workflow_output.__dict__)
            obj.sample_id = sample_id
            obj._parse()
            return obj

    def _parse(self):
        """Parse the variants annotated JSON file"""
        fmt = self.config.variants_annotated_json[self.workflow_id]
        self.path = Path(
            os.path.join(self.root, fmt.format(self.sample_id, self.sample_id))
        )
        if not self.path.is_file():
            logging.error(
                f"Small variant genome VCF missing: File {self.path} does not exist."
            )
            raise FileNotFoundError
        if self.path.suffix == ".gz":
            with gzip.open(self.path, "rt") as file:
                self.data = msgspec.json.decode(file.read())
        else:
            with open(self.path, "r") as file:
                self.data = msgspec.json.decode(file.read())
        self.variant_dict = self._get_variant_dict()

    def _get_variant_dict(self) -> dict[str, tuple[str | None, list[str]]]:
        """Create a dictionary mapping variant IDs to their type and overlapping genes from the annotated JSON data."""
        json_by_variant = {}
        for position in self.data.get("positions", []):
            # Get variant type
            variants = position.get("variants", [{}])
            variant_type = variants[0].get("variantType", None)  # if variants else None

            # Get overlapping genes
            genes = (
                sorted(
                    {
                        transcript["hgnc"]
                        for transcript in variants[0].get("transcripts", [])
                        if transcript.get("hgnc")
                    }
                )
                if variants
                else []
            )

            # Generate variant IDs for each alternate allele
            for alternate in position.get("altAlleles", []):
                variant_id = get_variant_id(
                    position["chromosome"],
                    position["position"],
                    position["refAllele"],
                    alternate,
                )
                json_by_variant[variant_id] = (variant_type, genes)
        return json_by_variant


def get_variant_id(
    chromosome: str, position: int, reference: str, alternate: str
) -> str:
    """Create the shared variant identifier used by the key variant sources."""
    return f"{chromosome}:{position}:{reference}>{alternate}"
