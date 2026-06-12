import logging
import os
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

    def workflow_id(self):
        """Return combined string for workflow type and version."""
        return f"{self.workflow_type}_{self.workflow_version}"


class SmallVariantGenomeVcf(WorkflowOutput):
    """Input class for small variant genome VCF files produced by different workflows.

    Attributes:
        path: Path to vcf (Path)
        sample_id: Sample identifier (str)
        vcf: Parsed VCF object (cyvcf2.VCF)
    """

    def __init__(self, config_yaml: str | Path, root_path: str | Path, sample_id: str):
        """Initialize SmallVariantGenomeVcf"""
        super().__init__(config_yaml, root_path)
        self.sample_id = sample_id

    @classmethod
    def create(cls, workflow_output: WorkflowOutput, sample_id: str):
        """Create SmallVariantGenomeVcf from existing WorkflowOutput"""
        obj = cls.__new__(cls)
        obj.__dict__.update(workflow_output.__dict__)
        obj.sample_id = sample_id
        return obj

    def parse(self):
        """Parse the VCF file"""
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
        return self.vcf


class TmbTrace(WorkflowOutput):
    """Input class for TMB trace files produced by different workflows.

    Attributes:
        path: Path to vcf (Path)
        rows: Parsed rows of the TMB trace file (polars.DataFrame)
        sample_id: Sample identifier (str)
    """

    def __init__(self, config_yaml: str | Path, root_path: str | Path, sample_id: str):
        """Initialize TmTraceTsv."""
        super().__init__(config_yaml, root_path)
        self.sample_id = sample_id

    @classmethod
    def create(cls, workflow_output: WorkflowOutput, sample_id: str):
        """Create TmbTraceTsv from existing WorkflowOutput."""
        obj = cls.__new__(cls)
        obj.__dict__.update(workflow_output.__dict__)
        obj.sample_id = sample_id
        return obj

    def parse(self):
        """Parse the TMB trace file"""
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
        return self.table


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

    @classmethod
    def create(cls, workflow_output: WorkflowOutput, sample_id: str):
        """Create VariantsAnnotatedJson from existing WorkflowOutput."""
        obj = cls.__new__(cls)
        obj.__dict__.update(workflow_output.__dict__)
        obj.sample_id = sample_id
        return obj

    def parse(self):
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
        with open(self.path, "r") as file:
            self.data = msgspec.json.decode(file.read())
        return self.data
