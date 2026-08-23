# ADA-PandA Quality Audit Log

## Purpose

This log records the independent quality review performed before revising the EJDR manuscript. The objective is to distinguish reproducibility from robustness, identify target leakage, detect inconsistent outputs, and prevent overclaiming.

## Source materials reviewed

- Submitted manuscript PDF: `EJDR-D-26-00174.pdf`
- Recovered notebook project: `ada-panda-mini`
- Notebooks:
  - `01_simulate_from_yaml.ipynb`
  - `02_compare_scenarios_UPDATED.ipynb`
  - `03_exposure_response.ipynb`
  - `04_benchmark_vs_literature.ipynb`
- Configuration: `configs/lit_anchored.yaml`
- Stored outputs:
  - `data/pk_cohort.csv`
  - `data/pk_cohort_with_detection.csv`
  - `reports/tlgs.parquet`
  - `reports/benchmarks.parquet`

## Finding 1: stored outputs are inconsistent

`data/pk_cohort.csv` contains 120 subjects but all records have `ADA_pos = 0`. It cannot be the source of the reported ADA-positive results.

`data/pk_cohort_with_detection.csv` contains ADA-positive records and produces the reported-style detection summaries, but its subject-level concentrations and conversion times do not match `pk_cohort.csv`. The stored project therefore contains outputs from different runs or notebook states.

Action: one authoritative pipeline and one manifest must replace the mixed outputs.

## Finding 2: target leakage through post-simulation scaling

Notebook 01 computes the Week 52 means and then rescales them to the configured targets:

```python
sf_neg = target_neg / mean_neg
sf_pos = target_pos / mean_pos
```

This makes the Week 52 means calibration anchors enforced by construction. They must not be presented as independent model validation. The primary analysis must remove this scaling. If retained, it must be labeled as a calibration/sensitivity scenario.

## Finding 3: manuscript/code mismatch

The manuscript identifies `scripts/run_pipeline.py` as the source-of-truth script, but the recovered project contains notebooks and no such script. The manuscript must either identify the notebooks as the source or add and validate a portable pipeline script.

## Finding 4: reproducibility is not robustness

The original notebooks use random conversion times, subject variability, residual error, and assay detection. A fixed seed reproduces one stochastic realization; it does not quantify sampling variability.

Correct terminology: “stochastic but reproducible for a fixed seed.”

## Finding 5: exact prevalence was forced

The original conversion function rounded anchor percentages to exact subject counts and sampled exactly those counts. For `n = 120`, the Week-52 true prevalence was therefore fixed at 49/120 = 40.83% in every replicate.

This is not a Monte Carlo result. It is a fixed-design constraint.

## Corrective analysis

A separate no-rescaling audit simulator was built outside the project first. It:

1. Reads the literature-anchored YAML configuration.
2. Converts cumulative ADA prevalence anchors into interval hazards.
3. Samples conversion independently for each subject.
4. Allows prevalence to vary across replicate cohorts.
5. Simulates one-compartment repeated-dose exposure.
6. Applies subject-level lognormal variability and residual error.
7. Applies ADA-related clearance after conversion without post-hoc target scaling.
8. Applies standard and PandA-like detection as visit-level observation processes.
9. Runs 1,000 independent seeds.
10. Reports replicate means, between-replicate standard deviations, Monte Carlo standard errors, and 2.5th/97.5th percentiles.

An intermediate draft omitted the interval before the first ADA anchor and produced implausibly low prevalence. That error was detected and corrected before accepting the results.

## Interim corrected results

Using 1,000 independent seeds and the current literature-anchored configuration:

| Quantity | Mean | MCSE | 2.5th–97.5th percentile |
|---|---:|---:|---:|
| True Week-52 ADA prevalence | 41.22% | 0.15 percentage points | 32.50–50.83% |
| Standard observed prevalence | 24.82% | 0.13 percentage points | 17.50–32.50% |
| PandA observed prevalence | 37.17% | 0.14 percentage points | 28.33–46.67% |
| Standard assay bias | −16.41 percentage points | 0.11 percentage points | −23.33 to −10.00 |
| PandA assay bias | −4.05 percentage points | 0.06 percentage points | −8.33 to −0.83 |
| True ADA-negative target attainment | 71.70% | 0.17 percentage points | 61.04–82.19% |
| Standard-observed ADA-negative attainment | 67.24% | 0.16 percentage points | 57.60–76.60% |
| PandA-observed ADA-negative attainment | 70.27% | 0.16 percentage points | 59.99–80.52% |

These are interim audit results, not final manuscript results. The simulator still requires code review, validation tests, and a fully portable implementation before publication use.

## Current interpretation

The corrected analysis supports an illustrative direction:

- standard drug-sensitive observation tends to under-detect ADA;
- PandA-like observation reduces, but does not eliminate, under-detection;
- observed ADA-negative target attainment can be lower than true ADA-negative target attainment.

It does not yet justify claims of a robust, product-specific, or clinically predictive effect.

## Required next steps

- Implement the corrected simulator as a portable project script.
- Add unit tests for anchor handling, conversion persistence, scaling removal, and assay observation.
- Decide whether the primary model uses interval hazards or a more explicit time-to-ADA model.
- Add replicate-level figures and uncertainty intervals.
- Add a separate calibration/sensitivity analysis only if scientifically justified.
- Pin the GitHub commit used for the manuscript.
- Rebuild tables and figures from one clean manifest.
- Rewrite the manuscript Results and Discussion conservatively.
- Draft the response letter only after the revised manuscript is complete.
