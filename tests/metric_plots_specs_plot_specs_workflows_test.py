"""Unit tests pinning production PLOT_SPECS content."""

from tsoppy.metric_plots.specs.plot_specs_workflows import PLOT_SPECS


def test_dragen_contamination_scatter_is_disabled():
    """LocalApp contamination scatters remain disabled for DRAGEN."""
    for spec_name in (
        "DNA_CONTAMINATION_P_VALUE_BY_RUN",
        "DNA_CONTAMINATION_P_VALUE_HIGHLIGHTED_RUN",
    ):
        assert PLOT_SPECS[spec_name]["dragen"] == {
            "plot": False,
            "index": 0,
        }


def test_dragen_contamination_score_bar_is_enabled():
    """DRAGEN uses contamination-score bar plot at index 20."""
    spec = PLOT_SPECS["DNA_CONTAMINATION_SCORE"]

    assert spec["dragen"] == {
        "plot": True,
        "index": 20,
    }

    assert spec["localapp"] == {
        "plot": False,
        "index": 0,
    }

    assert spec["plot_kind"] == "bar"

    assert spec["y_var"] == "DNA_CONTAMINATION_SCORE"


def test_dragen_run_level_plots_are_enabled_in_order():
    """Five shared run-level plots remain DRAGEN indices 1-5."""
    expected = {
        "CLUSTERS_PASSING_FILTER": 1,
        "ESTIMATED_YIELD": 2,
        "PCT_PF_READS": 3,
        "PCT_Q30_R1": 4,
        "PCT_Q30_R2": 5,
    }

    for (
        spec_name,
        expected_index,
    ) in expected.items():
        assert PLOT_SPECS[spec_name]["dragen"] == {
            "plot": True,
            "index": expected_index,
        }


def test_all_enabled_plot_indices_are_positive():
    """Enabled plots always have positive indices."""
    for workflow in (
        "dragen",
        "localapp",
    ):
        for spec in PLOT_SPECS.values():
            if spec[workflow]["plot"]:
                assert spec[workflow]["index"] > 0


def test_all_disabled_plot_indices_are_zero():
    """Disabled plots consistently use index zero."""
    for workflow in (
        "dragen",
        "localapp",
    ):
        for spec in PLOT_SPECS.values():
            if not spec[workflow]["plot"]:
                assert spec[workflow]["index"] == 0
