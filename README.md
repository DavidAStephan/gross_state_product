# Quarterly Gross State Product for Australian states and territories

The ABS publishes GSP annually and State Final Demand (SFD) quarterly. This project estimates
seasonally adjusted **quarterly GSP for all eight states and territories**, in current prices and
chain volume measures, by estimating the gap between GSP and SFD with two methods:

- **Method A** – temporal disaggregation of the annual gap on quarterly indicators, then
  movement-preserving two-way reconciliation (Di Fonzo–Marini).
- **Method B** – a joint linear Gaussian state-space model (KFAS) with 68% and 90% intervals.
- **Method C** (alternate) – the mixed-frequency VAR with stochastic volatility and Dirichlet-Laplace
  shrinkage of Koop, McIntyre, Mitchell and Poon (2020), with each posterior draw reconciled to the
  constraints. Less accurate than A and B on Australian data; see the report.

Both satisfy, to a relative tolerance of 1e-6 (in practice about 1e-12): each state's four quarters
(Sep–Jun) sum to annual GSP, and the states sum each quarter to the national target.

Start with `output/report.html`: method in plain language, charts, validation tables,
recommendation, and every judgement call.

## How to rerun when new ABS data are released

```sh
Rscript -e 'renv::restore()'   # first time only: installs the recorded package versions
Rscript run_all.R              # downloads the latest ABS data, runs 01-06, renders the report
```

That is all that is normally needed: nothing in the code is tied to a particular quarter or year.
The run takes about 50 minutes (Method C and its real-time test are the slowest parts). Then check, in this order:

1. **The console output of `01_data.R`** (or `output/check_sa_sums.csv` and
   `output/check_constraint_consistency.csv`): the sample dates, the number of unbenchmarked
   quarters, and how far seasonally adjusted quarters are from annual totals.
2. **`output/report.html`**: tests passed, divergence flags for the latest quarters.

What changes with each release:

| Release | Effect |
|---|---|
| Quarterly national accounts (March, June, September, December) | One more quarter. Between November and the next September-quarter accounts there are 1–4 unbenchmarked quarters; these have only the cross-sectional constraint. |
| Annual state accounts (November) | A new benchmark year. The previous nowcast year becomes benchmarked, industry shares roll forward one year, and model selection and the real-time test move forward one year. |
| Monthly trade and labour force | Used only up to the last complete quarter of the national accounts. |

Run the pipeline after the **quarterly national accounts**, since that release sets the last
quarter. Running just after the annual state accounts in November but before the December release
of the September-quarter national accounts works too (it simply has no unbenchmarked quarters).

### Options (environment variables)

| Variable | Default | Meaning |
|---|---|---|
| `GSP_REFRESH` | `true` | `false` reuses the spreadsheets already in `data/raw/` instead of downloading. |
| `GSP_USE_BU` | `false` | `true` lets Method A's model selection use the bottom-up industry indicator sets (tested in `04b_bottom_up_test.R`; they did not improve accuracy). |
| `GSP_TARGET` | `denton` | How national GDP is adjusted to form the cross-sectional target. `denton` interpolates the adjustment ratio smoothly; `prorata` applies one ratio per financial year (a step each July). Annual totals are identical. |

### If something fails

The scripts stop with a message starting `DATA PROBLEM:` that names the table or file needed. They
never substitute placeholder data.

- **ABS download fails**: download the named time series spreadsheets from the ABS website into
  `data/raw/` and rerun with `GSP_REFRESH=false`.
- **"expected exactly one series ... but found 0"**: the ABS renamed a series or table. Each series is
  selected in `01_data.R` by a table pattern and a series-name pattern; fix the pattern there. The
  current matches are in `output/data_dictionary.csv`.
- **Queensland Treasury workbook**: publication is suspended and the file is no longer on the live
  QGSO site. The script tries the September quarter 2020 issue (original address, then the Internet
  Archive) and falls back to the March quarter 2019 issue on `archive.qgso.qld.gov.au`. To use a
  specific issue, put `qld-state-accounts-tables-202009.xlsx` or `...-201903.xlsx` in
  `data/raw/qld/`. A workbook already there is used as-is.
- **A test fails in `06_validate.R`**: the run stops. The constraint checks in `04_method_a.R` and
  `05_method_b.R` would normally have stopped it earlier; look there first.

## Files

| File | Purpose |
|---|---|
| `01_data.R` | Downloads ABS data (`readabs`), builds tidy inputs, writes the data dictionary, checks seasonally adjusted sums and the consistency of the two constraint sets, builds the national target, fetches the Queensland Treasury workbook. |
| `02_indicators.R` | Production-side indicator (national industry value added × lagged state industry shares, chain-linked), indicator panel, annual gap. |
| `03_baselines.R` | Pro-rata and Denton on SFD, scaled to the national total. |
| `04_method_a.R` | Out-of-sample errors for every candidate model, model selection per state, preliminary estimates, reconciliation, constraint check. The slowest step (about 5 minutes). |
| `04b_bottom_up_test.R` | Tests the bottom-up industry indicator in Methods A and B against the lagged-share indicator (pseudo real-time and Queensland Treasury). |
| `05_method_b.R` | State-space model, smoothed estimates and intervals, constraint check. |
| `05b_method_c.R` | MF-VAR-SV (Gibbs sampler, about 15 minutes; needs the `GIGrvg` package). Optional: if it has not been run, `06_validate.R` leaves Method C out. |
| `06_validate.R` | Writes `output/quarterly_gsp.csv`, runs the `testthat` tests, pseudo real-time test, Queensland cross-check, A-versus-B comparison, sanity checks. |
| `report.qmd` | Quarto report (reads the saved results; does no estimation). |
| `07_dashboard_export.R` | Optional: small JSON files in `output/dashboard/` for the results dashboard. |
| `run_all.R` | Runs everything in order. |
| `R/utils.R` | Constants, date helpers. |
| `R/data_helpers.R` | Download, series selection, seasonal adjustment of trade, Queensland reader. |
| `R/constraints.R` | National target, pseudo real-time vintages, constraint checks. |
| `R/indicators.R` | Production indicator. |
| `R/baselines.R`, `R/method_a.R`, `R/method_b.R`, `R/method_c.R` | The estimators. Each takes the same inputs list, so they can be run on any vintage. |
| `R/validate.R` | Pseudo real-time loop and comparison helpers. |
| `tests/testthat/test-constraints.R` | Constraint and output tests. |
| `renv.lock` | Package versions. |

Scripts must be run from the project root and in order; each saves its results to
`data/processed/` for the next.

## Outputs (`output/`)

| File | Contents |
|---|---|
| `quarterly_gsp.csv` | Long format: `state`, `quarter` (first day of the quarter's last month, e.g. `2026-06-01` = June quarter 2026), `price_basis` (`current`, `cvm`), `method` (`method_a`, `method_b`, `method_c`, `baseline_prorata`, `baseline_denton`), `estimate` ($m, seasonally adjusted), `lower_68`, `upper_68`, `lower_90`, `upper_90` (Methods B and C), `benchmarked` (`FALSE` = nowcast quarter). |
| `report.html` | The report. |
| `data_dictionary.csv` | Every ABS series used: catalogue, table, series name, type, series ID, units, span. |
| `check_sa_sums.csv`, `check_constraint_consistency.csv` | The pre-estimation checks. |
| `method_a_*.csv`, `method_b_*.csv` | Model selection, coefficients, reconciliation adjustments, state-space parameters and uncertainty. |
| `validation_*.csv`, `sanity_*.csv` | Real-time errors, Queensland comparison, A-versus-B divergences, growth flags, gap shares. |

## Extending

- **Add an indicator**: add the series to the `q_state` block in `01_data.R` (named
  `<name>_current` and `<name>_cvm`), add `<name>` to the `select()` in `build_panel()` in
  `R/indicators.R`, then add a specification to `SPECS` in `R/method_a.R`. It will be included in
  model selection.
- **Change the indicators in Method B**: `ssm_data()` in `R/method_b.R` builds the two indicator
  matrices `Z1` and `Z2`.
- **Run any method on an earlier vintage**: `vintage_inputs(load_inputs(), last_fy)` followed by
  `method_a()`, `method_b()` or `run_baselines()`.
- Method C settings (lags, chain length, prior hyperparameters) are in `MFVAR_DEFAULT` in `R/method_c.R`; the header lists where it departs from the paper.

## Known limitations

See "Judgement calls and caveats" in the report. In brief: state trade indicators are goods only
and seasonally adjusted here rather than by the ABS; the within-year quarterly profile is much less
certain than the annual totals; NT and ACT estimates are volatile; Method B's intervals understate
uncertainty; the real-time test uses today's data vintage.
