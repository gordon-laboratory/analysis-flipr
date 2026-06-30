# CLAUDE.md — analysis-flipr

## Project overview
Analysis pipeline for **lifetime photometry (FLIPR)** data.

- `scripts_preprocessing/` — raw data preprocessing scripts (signal extraction, motion correction, normalization, etc.)
- `scripts_analysis/` — exploratory, visualization, summary, and statistics scripts
- `data/` — working data (not tracked by git; add to .gitignore)
- `data_example/` — small example datasets for testing/demo
- `output/` — saved results (figures, tables, stats); subdirectories mirror script names
- `prompts/` — AI prompt logs used during development

## GitHub
Organization: **gordon-laboratory**
Repo: `gordon-laboratory/analysis-flipr`

Git user for this repo is `gordon-laboratory`. Before pushing, ensure this is the
active `gh` account (`gh auth status`). If git picks up the wrong credential, run:
```bash
gh auth setup-git
git push origin master
```

## Data layout

```
data/sessions/<session>/data_flipr/   ← all preprocessed outputs
```

Session name is derived from the filename stem by stripping the trailing `_###` run number
(e.g. `2024_04_21_acz02_001.iFLiP2` → session `2024_04_21_acz02`).

All output CSVs have **no header row**. Column layout is documented in README.md.

## Functions

`functions/` contains shared Python modules ported from `source/iFLiP2_Code for Users/`.
Do not add ad-hoc MPET or lifetime calculations — always use these functions.

| Module | Exports |
|---|---|
| `read_iflip2.py` | `read_iflip2(filepath)` → dict with `data`, `header`, `LTTime`, `sampleTime`, `marks` |
| `calculate_mpet.py` | `calculate_mpet(data, header, spc_range, t0, ...)` → `(mpet, corrected_data)` |
| `fit_lifetime.py` | `fit_single_exp(lt_time, intensity, ...)`, `fit_double_exp(lt_time, intensity, ...)` |

## Preprocessing scripts

`scripts_preprocessing/01_tidy_and_preprocess.py` — single entry point for all preprocessing.
Takes a `.iFLiP2` file; requires `--spc_range` and `--t0`. See README for full argument list.

## Conventions
- Python scripts; prefer `pandas`, `numpy`, `matplotlib`/`seaborn` for core work
- Preprocessing scripts should be runnable standalone with clear CLI arguments
- Keep `data/` out of version control (add `.gitignore` entry); commit `data_example/` only

## Working with Claude
- Consult `CLAUDE_log.md` for session history before starting new work
- Update `CLAUDE_log.md` at the end of each session with a brief summary
