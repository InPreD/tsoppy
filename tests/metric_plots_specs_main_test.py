"""Unit tests pinning production PLOT_SPECS content."""

from pytest import mark

from tsoppy.metric_plots.specs.main import PLOT_SPECS


@mark.parametrize(
    "inputs, want",
    [
        # LocalApp contamination scatters remain disabled for DRAGEN.
        (
            ("DNA_CONTAMINATION_P_VALUE_BY_RUN", "dragen"),
            {
                "plot": False,
                "index": 0,
            },
        ),
        (
            ("DNA_CONTAMINATION_P_VALUE_HIGHLIGHTED_RUN", "dragen"),
            {
                "plot": False,
                "index": 0,
            },
        ),
        # DRAGEN uses contamination-score bar plot at index 20.
        (
            ("DNA_CONTAMINATION_SCORE", "dragen"),
            {
                "plot": True,
                "index": 20,
            },
        ),
        (
            ("DNA_CONTAMINATION_SCORE", "localapp"),
            {
                "plot": False,
                "index": 0,
            },
        ),
        (("DNA_CONTAMINATION_SCORE", "plot_kind"), "bar"),
        (("DNA_CONTAMINATION_SCORE", "y_var"), "DNA_CONTAMINATION_SCORE"),
        # Five shared run-level plots remain DRAGEN indices 1-5.
        (
            ("CLUSTERS_PASSING_FILTER", "dragen"),
            {
                "plot": True,
                "index": 1,
            },
        ),
        (
            ("ESTIMATED_YIELD", "dragen"),
            {
                "plot": True,
                "index": 2,
            },
        ),
        (
            ("PCT_PF_READS", "dragen"),
            {
                "plot": True,
                "index": 3,
            },
        ),
        (
            ("PCT_Q30_R1", "dragen"),
            {
                "plot": True,
                "index": 4,
            },
        ),
        (
            ("PCT_Q30_R2", "dragen"),
            {
                "plot": True,
                "index": 5,
            },
        ),
    ],
)
def test_plot_specs_specific(inputs, want):
    assert PLOT_SPECS[inputs[0]][inputs[1]] == want


def test_plot_specs_all():
    """Check that all plot specs have valid indices for workflows that enable plotting."""
    for workflow in (
        "dragen",
        "localapp",
    ):
        for spec in PLOT_SPECS.values():
            if spec[workflow]["plot"]:
                assert spec[workflow]["index"] > 0
            else:
                assert spec[workflow]["index"] == 0
