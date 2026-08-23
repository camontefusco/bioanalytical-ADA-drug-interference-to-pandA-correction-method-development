from pathlib import Path
import importlib.util

import numpy as np
import yaml


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("run_audit", ROOT / "scripts" / "run_audit.py")
run_audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(run_audit)
CFG = yaml.safe_load((ROOT / "configs" / "lit_anchored.yaml").read_text())


def test_fixed_seed_is_reproducible():
    first = run_audit.simulate(12345, CFG)
    second = run_audit.simulate(12345, CFG)
    assert first == second


def test_independent_seeds_produce_variable_prevalence():
    values = [run_audit.simulate(seed, CFG)["true_prev"] for seed in range(100, 120)]
    assert len(set(values)) > 1
    assert np.std(values, ddof=1) > 0


def test_observed_positives_cannot_exceed_true_positives():
    for seed in range(10, 20):
        result = run_audit.simulate(seed, CFG)
        assert result["standard_prev"] <= result["true_prev"]
        assert result["panda_prev"] <= result["true_prev"]


def test_replicate_summary_reports_mcse():
    rows = [run_audit.simulate(seed, CFG) for seed in range(20)]
    summary = run_audit.summarize(
        __import__("pandas").DataFrame(rows)
    )
    assert summary["true_prev"]["mcse"] > 0
    assert summary["standard_bias_pp"]["q025"] <= summary["standard_bias_pp"]["q975"]


def test_expected_observation_estimands_are_bounded():
    result = run_audit.simulate(12345, CFG)
    assert 0 <= result["standard_expected_prev"] <= result["true_prev"]
    assert 0 <= result["panda_expected_prev"] <= result["true_prev"]


def test_primary_run_has_no_target_rescaling():
    result = run_audit.simulate(12345, CFG)
    configured_target = CFG["trough_targets_week52"]["ADA_negative_ng_per_mL"]
    # The generated mean is stochastic and must not be forced to equal the target.
    assert not np.isclose(
        result["true_negative_mean_trough"],
        configured_target,
        rtol=0,
        atol=1e-9,
    )
