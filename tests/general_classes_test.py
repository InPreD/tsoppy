from os import path
from contextlib import nullcontext
from pytest import mark, raises

from tsoppy.general.classes import (
    WorkflowOutput,
    SmallVariantGenomeVcf
)

# Define path to test data - cannot be absolute due to different paths locally and in CI
test_data_dir = "tests/test_data/general_classes"


@mark.parametrize(
    "inputs, exception, want",
    [
        (
            (
                'config.yaml',
                path.join(test_data_dir, "dragen/standard")
            ),
            nullcontext(),
            'dragen_2.6.2.4'
        ),
        (
            (
                'config.yaml',
                path.join(test_data_dir, "localapp/standard")
            ),
            nullcontext(),
            'localapp_ruo-2.2.0.12'
        ),
        (
            (
                'config.yaml',
                path.join(test_data_dir, "localapp/non-existent")
            ),
            raises(FileNotFoundError),
            ''
        )
    ]
)
def test_workflowoutput_init(inputs, exception, want):
    with exception:
        got = WorkflowOutput(inputs[0], inputs[1])
        assert got.workflow_id() == want
