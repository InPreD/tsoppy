
from file_parser import Parse_section_tsv
import logging
import os
from pathlib import Path

import cyvcf2
import msgspec
import polars


# Use logger that was set up in CLI
logger = logging.getLogger(__name__)

class MetricsOutputTsv(WorkflowOutput):
    """Input class for MetricsOutput.tsv files produced by different workflows.

    Attributes:
        path: Path to MetricsOutput.tsv
        headers: Header lines parsed from the file
        sections: Parsed section DataFrames
        table: Parsed section dictionary
    """

    def __init__(self, config_yaml: str | Path, root_path: str | Path):
        """Initialize MetricsOutputTsv."""
        super().__init__(config_yaml, root_path)
        self._parse()

    @classmethod
    def create(cls, workflow_output: WorkflowOutput):
        """Create MetricsOutputTsv from existing WorkflowOutput."""
        obj = cls.__new__(cls)
        obj.__dict__.update(workflow_output.__dict__)
        obj._parse()
        return obj

    def _parse(self):
        """Parse the MetricsOutput.tsv file."""
        info_src = list(self.config.metrics_output_tsv.values())

        if len(set(info_src)) != 1:
            raise ValueError(
                f"Got {info_src} but need exactly one MetricsOutput.tsv path"
            )

        self.path = Path(os.path.join(self.root, info_src[0]))

        if not self.path.is_file():
            logging.error(
                f"MetricsOutput.tsv missing: File {self.path} does not exist."
            )
            raise FileNotFoundError

        self.headers, self.sections = Parse_section_tsv(
            str(self.path),
            ["Header"]
        )

        self.table = self.sections

    def get_section(self, section_name: str) -> polars.DataFrame:
        """Return one parsed MetricsOutput.tsv section."""
        if section_name not in self.sections:
            raise KeyError(
                f"Section '{section_name}' not found. "
                f"Available sections: {list(self.sections.keys())}"
            )
        return self.sections[section_name]

    def section_names(self) -> list[str]:
        """Return available MetricsOutput.tsv section names."""
        return list(self.sections.keys())
    

workflow = WorkflowOutput(
    config_yaml="workflow_config.yaml",
    root_path="/path/to/workflow/output" 
)

metrics = MetricsOutputTsv.create(workflow)

print(metrics.path)
print(metrics.workflow_type)
print(metrics.workflow_version)
print(metrics.section_names())

dna_qc = metrics.get_section("DNA Library QC Metrics")
print(dna_qc)