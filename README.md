| Methods | Code |
|---|---|
| Merge frailty, PGSs and ancestry PCs | Read and join three ID-keyed CSV inputs |
| Restrict ELSA participants to age 65+ | `min_age = 65` |
| 26 frailty-relevant PGSs | `n_pgs_expected = 26` |
| Square-root transformed FI | `sqrt_frailty = TRUE` |
| Standardised FI | `standardise_frailty = TRUE` |
| Single-PGS model with age, sex and 10 PCs | Null/full ordinary least squares and incremental R² |
| Joint multi-PGS prediction | Elastic net (`alpha = 0.5`) |
| Repeated tenfold cross-validation, ten repeats | `caret::trainControl(repeatedcv, 10, 10)` |
| Ten lambda candidates chosen by RMSE | `caret::train` over 10 candidate lambda values |
| Bootstrap coefficient stability | 1,000 resamples at CV-selected lambda |
| Nonlinear frailty progression | `gamm4` with `s(age)` and `s(age, by = PGS)` |
| Account for repeated observations | Participant-specific random intercept |
| At least three observations per participant | `minimum_long_visits = 3` |

### Inputs

- `frailty_long.csv`: `id, age, sex, wave, frailty` (one record per participant per wave)
- `pgs.csv`: `id, PGS_1, ..., PGS_26` (one record per person)
- `pcs.csv`: `id, PC_1, ..., PC_10` (one record per person)

Column names, PGS selection and preprocessing rules are configurable at the start of `analysis.R`. The script uses earliest eligible age as the baseline frailty observation, and complete cases for regression inputs. It also derives the FI scaling constants from baseline observations. For longitudinal trajectories the same baseline-derived constants are used for all waves.

### Run

Install the R packages `dplyr`, `caret`, `glmnet`, `gamm4`, and `mgcv`, then:

```bash
MPS_analysis.R
```

Outputs including single-PGS incremental R², CV tuning, final elastic-net coefficients, bootstrap selection frequencies/intervals, and longitudinal interaction tests.

### Validation and interpretation caveats

- **LBC1936 is an independent cohort.** 
- LBC1936 analysis using four PCs and 1,000 bootstraps, 
- The single-PGS R² is the difference in *in-sample* R² between full and covariate-only regression, not an external predictive R².
- Mixed model interaction smooth p-values are approximate. GAMM analyses of multiple scores receive BH/FDR correction.
- The code currently requires complete cases. Longitudinal models adjust for sex and ancestry PCs and include separate age/PGS smooths.
