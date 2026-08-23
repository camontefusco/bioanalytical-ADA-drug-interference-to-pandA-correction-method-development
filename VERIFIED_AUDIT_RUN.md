# Verified Monte Carlo Audit Run

Date: 2026-08-23

## Repository verification

The committed repository version was downloaded and executed locally:

- `scripts/run_audit.py`
- `tests/test_run_audit.py`
- `configs/lit_anchored.yaml`

Test result: **5 passed**.

Command:

```bash
python3 -m pytest -q tests/test_run_audit.py
python3 scripts/run_audit.py --config configs/lit_anchored.yaml --replicates 1000 --seed-start 100000
```

## Verified 1,000-replicate results

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
| True ADA-negative mean trough | 10,247.5 ng/mL | 11.3 ng/mL | 9,580.7–11,010.1 |
| True ADA-positive mean trough | 8,784.3 ng/mL | 15.2 ng/mL | 7,871.6–9,743.8 |

## New conclusion

After removal of post-hoc target scaling, the model does not reproduce the submitted ADA-positive Week-52 trough anchor of 4,855 ng/mL. The submitted exposure separation therefore depended on the prior scaling/calibration workflow.

Before manuscript revision, the model must either:

1. use a transparent, pre-specified PK calibration that produces the intended exposure separation without post-hoc outcome scaling; or
2. present the current analysis as a generic illustrative scenario and remove claims that it reproduces the adalimumab Week-52 exposure anchors.

The second option is statistically safer unless a defensible calibration procedure and independent validation target are established.
