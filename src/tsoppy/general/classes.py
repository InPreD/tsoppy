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


class VariantsAnnotatedJson(WorkflowOutput):
    """Input class for annotated JSON files produced by different workflows.

    Attributes:
        path: Path to vcf (Path)
        data: Parsed JSON data (dict)
        sample_id: Sample identifier (str)
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
        else:
            with open(self.path, "r") as file:
                self.data = msgspec.json.decode(file.read())
