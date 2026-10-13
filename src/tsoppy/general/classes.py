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

    def sample_meta(self, sample_id: str) -> dict[str:str]:
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
        path: Path to vcf (Path)
        sample_id: Sample identifier (str)
        vcf: Parsed VCF object (cyvcf2.VCF)
    """

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


class TmbTraceTsv(WorkflowOutput):
    """Input class for TMB trace files produced by different workflows.

    Attributes:
        path: Path to vcf (Path)
        table: Parsed rows of the TMB trace file (polars.DataFrame)
        sample_id: Sample identifier (str)
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


class VariantsAnnotatedJson(WorkflowOutput):
    """Input class for annotated JSON files produced by different workflows.

    Attributes:
        path: Path to vcf (Path)
        data: Parsed JSON data (dict)
        sample_id: Sample identifier (str)
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
