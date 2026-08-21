import gzip
import logging
import os
import re
from pathlib import Path
import cyvcf2
import msgspec
import polars

from tsoppy.general.file_parser import Parse_section_tsv

# Use logger that was set up in CLI
logger = logging.getLogger(__name__)


class WorkflowConfig(msgspec.Struct):
    """Config class for workflow output file path format strings"""

    metrics_output_tsv: dict[str, str]
    small_variant_genome_vcf: dict[str, str]
    tmb_trace_tsv: dict[str, str]
    variants_annotated_json: dict[str, str]

    def __eq__(self, other):
        if not isinstance(other, WorkflowConfig):
            return False
        if self.metrics_output_tsv != other.metrics_output_tsv:
            return False
        if self.small_variant_genome_vcf != other.small_variant_genome_vcf:
            return False
        if self.tmb_trace_tsv != other.tmb_trace_tsv:
            return False
        return self.variants_annotated_json == other.variants_annotated_json


class WorkflowOutput:
    """Base class for outputs produced by different workflows (e.g. dragen/localapp).

    Attributes:
        config: Configuration (WorkflowConfig)
        root: Root path (Path)
        workflow_type: Detected workflow type (str)
        workflow_version: Detected workflow version (str)
    """

    def __init__(self, config_yaml: str | Path, root_path: str | Path):
        """Initialize WorkflowOutput."""
        self.root = Path(root_path)
        with open(config_yaml, "r") as yaml_file:
            self.config = msgspec.yaml.decode(yaml_file.read(), type=WorkflowConfig)

        self._detect_type_and_version()

    def _detect_type_and_version(self):
        """Detect which workflow type and version is present based on information in MetricsOutput.tsv."""

        # Get all values for MetricsOutput.tsv paths and check if they are the same
        info_src = list(self.config.metrics_output_tsv.values())
        if len(set(info_src)) != 1:
            raise ValueError(
                f"Got {info_src} but need exactly one file to detect workflow id"
            )

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

    def __eq__(self, other):
        if not isinstance(other, WorkflowOutput):
            return False
        if self.config != other.config:
            return False
        if self.root != other.root:
            return False
        if self.workflow_type != other.workflow_type:
            return False
        return self.workflow_version == other.workflow_version

    def workflow_id(self):
        """Return combined string for workflow type and version."""
        return f"{self.workflow_type}_{self.workflow_version}"


class SmallVariantGenomeVcf(WorkflowOutput):
    """Input class for small variant genome VCF files produced by different workflows.

    Attributes:
        header_dict: Dict containing parsed VCF object header information accessible by keys (dict)
        header_rex: Regex to parse top level key value pairs in VCF object header (str)
        header_rex: Regex to parse sub level key value pairs in VCF object header (str)
        path: Path to vcf (Path)
        sample_id: Sample identifier (str)
        vcf: Parsed VCF object (cyvcf2.VCF)
    """

    header_dict = {}
    header_rex = r"##(?P<key>\w+)=(?P<value>.+)"
    header_subrex = r"(\w+)=([\w\d]+|\".+\")?"

    def __init__(self, config_yaml: str | Path, root_path: str | Path, sample_id: str):
        """Initialize SmallVariantGenomeVcf"""
        super().__init__(config_yaml, root_path)
        self.sample_id = sample_id
        self._parse()

    @classmethod
    def create(cls, workflow_output: WorkflowOutput, sample_id: str):
        """Create SmallVariantGenomeVcf from existing WorkflowOutput"""
        obj = cls.__new__(cls)
        obj.__dict__.update(workflow_output.__dict__)
        obj.sample_id = sample_id
        obj._parse()
        return obj

    @staticmethod
    def _variant_id(
        chromosome: str, position: int, reference: str, alternate: str
    ) -> str:
        """Create the shared variant identifier used by the key variant sources."""
        return f"{chromosome}:{position}:{reference}:{alternate}"

    @classmethod
    def merge(
        cls,
        workflow_output: WorkflowOutput,
        sample_id: str,
        output_path: str | Path | None = None,
    ) -> None | polars.DataFrame:
        """Merge VCF enriched with TMB and Nirvana annotations and either save to file or return polars dataframe."""

        vcf_obj = cls.create(workflow_output, sample_id)
        tmb_obj = TmbTraceTsv.create(workflow_output, sample_id)
        json_obj = VariantsAnnotatedJson.create(workflow_output, sample_id)

        if output_path:
            cls._write_vcf(vcf_obj, tmb_obj, json_obj, output_path)
        else:
            return cls._to_dataframe(vcf_obj, tmb_obj=tmb_obj, json_obj=json_obj)

    def _write_vcf(
        self,
        tmb_obj: TmbTraceTsv,
        json_obj: VariantsAnnotatedJson,
        output_path: str | Path,
    ):
        """Write a merged VCF with TMB and Nirvana annotations."""

        out_path = Path(output_path)
        tmb_by_variant = tmb_obj.variant_dict
        tmb_df = tmb_obj.class_table

        json_by_variant = json_obj.variant_dict
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
                "Description": "Variant classification provided by initial Illumina pipeline analysis: 'Germline_DB'', 'Germline_Proxi', 'Somatic', 'Blacklist', 'VCF_filtered'",
            },
        )

        for info in info_headers:
            self.vcf.add_info_to_header(info)

        writer = cyvcf2.Writer(out_path, self.vcf)
        for variant in self.vcf:
            if not variant.ALT or variant.ALT == ["<NON_REF>"]:
                continue
            variant_id = SmallVariantGenomeVcf._variant_id(
                variant.CHROM, variant.POS, variant.REF, variant.ALT[0]
            )
            variant.INFO["variant_ID"] = variant_id

            variant_type, genes = (
                json_by_variant[variant_id]
                if variant_id in json_by_variant
                else ("", [])
            )
            variant.INFO["variant_type"] = variant_type
            variant.INFO["overlapping_genes"] = ",".join(genes)
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

        print(f"Merged VCF file saved to '{out_path}'")

    def _to_dataframe(
        self, tmb_obj: TmbTraceTsv, json_obj: VariantsAnnotatedJson
    ) -> polars.DataFrame:
        """Return merged variant data as a Polars dataframe."""

        tmb_df = tmb_obj.class_table
        json_by_variant = json_obj.variant_dict
        var_stats = ["AD", "DP", "AF" if self.workflow_type == "dragen" else "VF"]
        cols = (
            ["variant_ID", "variant_type", "overlapping_genes"] + var_stats + ["filter"]
        )
        var_dict = {c: [] for c in cols}

        for variant in self.vcf:
            if not variant.ALT or variant.ALT == ["<NON_REF>"]:
                continue
            variant_id = self._variant_id(
                variant.CHROM, variant.POS, variant.REF, variant.ALT[0]
            )
            variant_type, genes = (
                json_by_variant[variant_id]
                if variant_id in json_by_variant
                else ("", [])
            )
            var_dict["variant_ID"].append(variant_id)
            var_dict["variant_type"].append(variant_type)
            var_dict["overlapping_genes"].append(genes)
            var_dict["filter"].append(variant.FILTER)
            for stat in var_stats:
                var_dict[stat].append(variant.format(stat)[0])

        df = polars.DataFrame(data=var_dict)

        # Check for single-element list columns
        explode_cols = [
            name
            for name, dtype in df.schema.items()
            if isinstance(dtype, polars.Array) and dtype.shape == (1,)
        ]
        if explode_cols:
            df = df.explode(explode_cols, empty_as_null=True)

        joined_df = df.join(
            tmb_df,
            left_on="variant_ID",
            right_on="variant_ID",
            how="left",
            maintain_order="left",
        )

        col_order = [
            "variant_ID",
            "variant_type",
            "Illumina_variant_class",
            "overlapping_genes",
        ] + var_stats

        if self.workflow_type == "dragen":
            full_df = (
                joined_df.with_columns(
                    polars.col("filter")
                    .fill_null(polars.col("Class"))
                    .replace({"excluded_regions": "Blacklist"}),
                )
                .drop("Class")
                .with_columns(
                    polars.when(polars.col("filter").str.contains(r"^[a-z]"))
                    .then(polars.lit("VCF_filtered"))
                    .otherwise(polars.col("filter"))
                    .alias("Illumina_variant_class")
                )
                .drop("filter")
                .select(col_order)
            )
        else:
            full_df = (
                joined_df.with_columns(
                    polars.when(polars.col("filter") == "Blacklist")
                    .then(polars.lit("Blacklist"))
                    .when(polars.col("filter").is_not_null())
                    .then(polars.lit("VCF_filtered"))
                    .otherwise(polars.col("filter"))
                    .alias("filter")
                )
                .with_columns(polars.col("filter").fill_null(polars.col("Class")))
                .rename({"filter": "Illumina_variant_class"})
                .drop("Class")
                .select(col_order)
            )

        return full_df

    def _parse(self):
        """Parse the small variant genome VCF file"""
        fmt = self.config.small_variant_genome_vcf[self.workflow_id()]
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
        """Parse vcf header into dict"""
        for match in re.finditer(self.header_rex, self.vcf.raw_header):
            # Use named capture groups to split line into key and value
            key = match.groupdict()["key"]
            value = match.groupdict()["value"]

            # Check if value containes additional key value pairs
            if not bool(re.search(self.header_subrex, value)):
                self.header_dict[key] = value
            else:
                # Parse key value pairs into dict
                item = dict()
                for match2 in re.finditer(self.header_subrex, value):
                    g = match2.groups()
                    item[g[0]] = g[1]

                # Use ID value as subkey if it exists
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

    def __init__(self, config_yaml: str | Path, root_path: str | Path, sample_id: str):
        """Initialize TmTraceTsv."""
        super().__init__(config_yaml, root_path)
        self.sample_id = sample_id
        self._parse()

    @classmethod
    def create(cls, workflow_output: WorkflowOutput, sample_id: str):
        """Create TmbTraceTsv from existing WorkflowOutput."""
        obj = cls.__new__(cls)
        obj.__dict__.update(workflow_output.__dict__)
        obj.sample_id = sample_id
        obj._parse()
        return obj

    def _parse(self):
        """Parse the TMB trace tsv."""
        fmt = self.config.tmb_trace_tsv[self.workflow_id()]
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
            SmallVariantGenomeVcf._variant_id(
                row["Chromosome"], row["Position"], row["RefCall"], row["AltCall"]
            ): row
            for row in self.table.iter_rows(named=True)
        }

    def _get_variant_class_table(self):
        tmb_df = self.table.with_columns(
            polars.concat_str(
                ["Chromosome", "Position", "RefCall", "AltCall"], separator=":"
            ).alias("variant_ID")
        )
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

    def __init__(self, config_yaml: str | Path, root_path: str | Path, sample_id: str):
        """Initialize VariantsAnnotatedJson."""
        super().__init__(config_yaml, root_path)
        self.sample_id = sample_id
        self._parse()

    @classmethod
    def create(cls, workflow_output: WorkflowOutput, sample_id: str):
        """Create VariantsAnnotatedJson from existing WorkflowOutput."""
        obj = cls.__new__(cls)
        obj.__dict__.update(workflow_output.__dict__)
        obj.sample_id = sample_id
        obj._parse()
        return obj

    def _parse(self):
        """Parse the variants annotated JSON file"""
        fmt = self.config.variants_annotated_json[self.workflow_id()]
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
                self.variant_dict = self._get_variant_dict()
        else:
            with open(self.path, "r") as file:
                self.data = msgspec.json.decode(file.read())
                self.variant_dict = self._get_variant_dict()

    def _get_variant_dict(self) -> dict[str, tuple[str | None, list[str]]]:
        json_by_variant = {}
        for position in self.data.get("positions", []):
            variants = position.get("variants", [{}])
            variant_type = variants[0].get("variantType") if variants else None
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
            for alternate in position.get("altAlleles", []):
                variant_id = SmallVariantGenomeVcf._variant_id(
                    position["chromosome"],
                    position["position"],
                    position["refAllele"],
                    alternate,
                )
                json_by_variant[variant_id] = (variant_type, genes)
        return json_by_variant
