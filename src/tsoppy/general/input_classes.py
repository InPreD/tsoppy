import os
from pathlib import Path
from typing import Dict, Optional

import cyvcf2
import msgspec
import polars


class BaseInput:
    """Base class for inputs produced by different workflows (e.g. dragen/localapp).

    Subclasses should define `default_subpath_formats` mapping workflow names to
    subpath format strings that accept `(sample, sample)` for formatting.

    Attributes:
        paths: Mapping of workflow names to resolved Path objects (Dict[str, Path])
        root: Root path (Path)
        sample: Sample name (str)
        subpath_formats: Mapping of workflow names to subpath format strings (Dict[str, str])
        type: Detected workflow type (str)
    """

    subpath_formats: Dict[str, str] = {}
    type: Optional[str] = None
    path: Optional[Path] = None

    def __init__(self, sample: str, root_path: str | Path, subpath_formats: Optional[Dict[str, str]] = None):
        """Initialize the BaseInput."""
        self.sample = sample
        self.root = Path(root_path)
        if subpath_formats:
            self.subpath_formats.update(subpath_formats)

        self._resolve_paths()
        self._detect_type()

    def _resolve_paths(self):
        """Resolve the paths for each workflow type based on the provided subpath formats."""
        out: Dict[str, Path] = {}
        for name, fmt in self.subpath_formats.items():
            out[name] = Path(os.path.join(
                self.root, fmt.format(self.sample, self.sample)))
        self.paths = out

    def _detect_type(self):
        """Detect which workflow type is present based on the existence of the resolved paths."""
        found = [name for name, path in self.paths.items() if path.is_file()]
        if len(found) > 1:
            raise ValueError(
                f"Multiple workflow files found for sample {self.sample}: {found}")
        if not found:
            raise FileNotFoundError(
                f"No workflow file found for sample {self.sample}. Searched: {self.paths}")
        self.type = found[0]
        self.path = self.paths[self.type]


class Vcf(BaseInput):
    """Input class for VCF files produced by different workflows.

    Attributes:
        default_subpath_formats: Mapping of workflow names to subpath format strings (Dict[str, str])
        vcf: Parsed VCF object (cyvcf2.VCF)
    """

    default_subpath_formats = {
        "dragen": "Logs_Intermediates/DnaDragenCaller/{}/{}.hard-filtered.gvcf.gz",
        "localapp": "Logs_Intermediates/VariantMatching/{}/{}_MergedSmallVariants.genome.vcf",
    }

    def __init__(self, sample: str, root_path: str | Path, subpath_formats: Optional[Dict[str, str]] = None):
        if subpath_formats:
            super().__init__(sample, root_path, subpath_formats)
        else:
            super().__init__(sample, root_path, self.default_subpath_formats)

    def parse(self):
        """Parse the VCF file"""
        self.vcf = cyvcf2.VCF(self.path)
        return self.vcf


class TmbTrace(BaseInput):
    """Input class for TMB trace files produced by different workflows.

    Attributes:
        rows: Parsed rows of the TMB trace file (polars.DataFrame)
    """

    default_subpath_formats = {
        "dragen": "Logs_Intermediates/Tmb/{}/{}.tmb.trace.tsv",
        "localapp": "Logs_Intermediates/Tmb/{}/{}_TMB_Trace.tsv",
    }

    def __init__(self, sample: str, root_path: str | Path, subpath_formats: Optional[Dict[str, str]] = None):
        if subpath_formats:
            super().__init__(sample, root_path, subpath_formats)
        else:
            super().__init__(sample, root_path, self.default_subpath_formats)

    def parse(self):
        """Parse the TMB trace file"""
        self.table = polars.read_csv(
            self.path, separator="\t")
        return self.table


class AnnotatedJson(BaseInput):
    """Input class for annotated JSON files produced by different workflows.

    Attributes:
        data: Parsed JSON data (dict)
    """

    default_subpath_formats = {
        "dragen": "Logs_Intermediates/Annotation/{}/{}_DNAVariants_Annotated.json",
        "localapp": "Logs_Intermediates/Annotation/{}/{}_SmallVariants_Annotated.json.gz",
    }

    def __init__(self, sample: str, root_path: str | Path, subpath_formats: Optional[Dict[str, str]] = None):
        if subpath_formats:
            super().__init__(sample, root_path, subpath_formats)
        else:
            super().__init__(sample, root_path, self.default_subpath_formats)

    def parse(self):
        """Parse the annotated JSON file"""
        with open(self.path, 'r') as file:
            self.data = msgspec.json.decode(file.read())
        return self.data
