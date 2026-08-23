from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def recovery(drug_ug_ml: np.ndarray, knots: list[list[float]]) -> np.ndarray:
    pairs = sorted((float(x), float(y)) for x, y in knots)
    x = np.array([p[0] for p in pairs])
    y = np.array([p[1] for p in pairs])
    d = np.clip(np.asarray(drug_ug_ml, float), x.min(), x.max())
    return np.interp(np.log10(d), np.log10(x), y)


def sample_conversion_weeks(
    rng: np.random.Generator,
    n: int,
    visits: np.ndarray,
    anchors: np.ndarray,
) -> np.ndarray:
    # Convert cumulative prevalence anchors to interval hazards.
    # Counts are sampled; they are not forced to equal the anchors.
    cumulative = np.maximum.accumulate(anchors)
    conversion = np.full(n, np.inf)
    previous = 0.0
    left_edges = np.concatenate(([0.0], visits[:-1]))
    for left, right, target in zip(left_edges, visits, cumulative):
        hazard = 0.0 if previous >= 1 else (target - previous) / (1.0 - previous)
        hazard = float(np.clip(hazard, 0.0, 1.0))
        at_risk = np.isinf(conversion)
        chosen = at_risk & (rng.random(n) < hazard)
        conversion[chosen] = left + rng.beta(2.0, 2.0, int(chosen.sum())) * (right - left)
        previous = float(target)
    return conversion


def simulate(seed: int, cfg: dict, n: int | None = None) -> dict:
    rng = np.random.default_rng(seed)
    n = int(n or cfg["cohort"]["n_subjects"])

    dose = float(cfg["dose_regimen"]["dose_mg"])
    tau = float(cfg["dose_regimen"]["interval_days"])
    bioavailability = float(cfg["pk_model"]["F"])
    volume = float(cfg["pk_model"]["V_L"])
    clearance = float(cfg["pk_model"]["CL_base_L_per_day"])
    ada_clearance_multiplier = float(cfg["pk_model"]["CL_multiplier_ADApos"])
    lloq = float(cfg["trough_targets_week52"]["lloq_ng_per_mL"])
    visits = np.asarray(cfg["ada_incidence"]["visits_weeks"], float)
    anchors = np.asarray(cfg["ada_incidence"]["pct_ADA_positive"], float) / 100.0

    times = np.arange(0.0, 365.0, 7.0)
    conversion = sample_conversion_weeks(rng, n, visits, anchors)
    true_ada = np.isfinite(conversion)
    concentrations = np.empty((n, len(times)))

    for i in range(n):
        base_clearance = clearance * (ada_clearance_multiplier if true_ada[i] else 1.0)
        k = base_clearance / volume
        subject = np.zeros(len(times))
        for dose_number in range(int(np.ceil(times.max() / tau)) + 2):
            relative_time = times - dose_number * tau
            mask = relative_time >= 0
            subject[mask] += (bioavailability * dose / volume) * np.exp(-k * relative_time[mask])

        if true_ada[i]:
            pre_conversion = times < conversion[i] * 7.0
            baseline = np.zeros(len(times))
            for dose_number in range(int(np.ceil(times.max() / tau)) + 2):
                relative_time = times - dose_number * tau
                mask = relative_time >= 0
                baseline[mask] += (bioavailability * dose / volume) * np.exp(
                    -(clearance / volume) * relative_time[mask]
                )
            subject[pre_conversion] = baseline[pre_conversion]

        concentrations[i] = subject * 1000.0

    iiv_sd = np.where(
        true_ada,
        float(cfg["cohort"]["iiv_lognormal_sd_ADApos"]),
        float(cfg["cohort"]["iiv_lognormal_sd_ADAneg"]),
    )
    concentrations *= np.exp(rng.normal(0, 1, (n, 1)) * iiv_sd[:, None])
    concentrations *= 1 + rng.normal(
        0, float(cfg["cohort"]["proportional_residual_sd"]), concentrations.shape
    )
    concentrations = np.clip(concentrations, 0, None)

    week52 = concentrations[:, -1]
    drug_ug_ml = week52 / 1000.0
    standard_probability = true_ada * recovery(
        drug_ug_ml, cfg["assay_recovery"]["standard"]["knot_points"]
    )
    panda_probability = true_ada * recovery(
        drug_ug_ml, cfg["assay_recovery"]["panda"]["knot_points"]
    )
    standard = (rng.random(n) < standard_probability).astype(int)
    panda = (rng.random(n) < panda_probability).astype(int)

    target = 0.8 * float(cfg["trough_targets_week52"]["ADA_negative_ng_per_mL"])
    return {
        "seed": seed,
        "n": n,
        "true_prev": float(true_ada.mean()),
        "standard_prev": float(standard.mean()),
        "panda_prev": float(panda.mean()),
        "standard_bias_pp": float((standard.mean() - true_ada.mean()) * 100),
        "panda_bias_pp": float((panda.mean() - true_ada.mean()) * 100),
        "true_negative_target_attainment": float((week52[~true_ada] >= target).mean()),
        "standard_observed_negative_target_attainment": float((week52[standard == 0] >= target).mean()),
        "panda_observed_negative_target_attainment": float((week52[panda == 0] >= target).mean()),
        "true_negative_mean_trough": float(week52[~true_ada].mean()),
        "true_positive_mean_trough": float(week52[true_ada].mean()),
        "lloq_fraction_week52": float((week52 < lloq).mean()),
    }


def summarize(results: pd.DataFrame) -> dict:
    summary = {}
    for column in results.columns:
        if column == "seed":
            continue
        values = results[column].astype(float)
        summary[column] = {
            "mean": float(values.mean()),
            "sd_between_replicates": float(values.std(ddof=1)),
            "mcse": float(values.std(ddof=1) / np.sqrt(len(values))),
            "q025": float(values.quantile(0.025)),
            "q975": float(values.quantile(0.975)),
        }
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the no-rescaling ADA-PandA Monte Carlo audit.")
    parser.add_argument("--config", default="configs/lit_anchored.yaml")
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed-start", type=int, default=100000)
    parser.add_argument("--subjects", type=int, default=None)
    parser.add_argument("--output-dir", default="reports/monte_carlo_audit")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    cfg = yaml.safe_load((root / args.config).read_text())
    output_dir = root / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    results = pd.DataFrame(
        [simulate(args.seed_start + i, cfg, args.subjects) for i in range(args.replicates)]
    )
    summary = summarize(results)
    results.to_csv(output_dir / "replicate_results.csv", index=False)
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    (output_dir / "run_manifest.json").write_text(
        json.dumps(
            {
                "config": args.config,
                "replicates": args.replicates,
                "seed_start": args.seed_start,
                "subjects": args.subjects or cfg["cohort"]["n_subjects"],
                "rescaling": False,
                "target_leakage_control": "Week-52 concentrations are not post-hoc scaled to target means.",
            },
            indent=2,
        )
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
