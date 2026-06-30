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

## Conventions
- Python scripts; prefer `pandas`, `numpy`, `matplotlib`/`seaborn` for core work
- Preprocessing scripts should be runnable standalone with clear CLI arguments
- Output filenames should include the script name and a timestamp or version tag
- Keep `data/` out of version control (add `.gitignore` entry); commit `data_example/` only

## Working with Claude
- Consult `CLAUDE_log.md` for session history before starting new work
- Update `CLAUDE_log.md` at the end of each session with a brief summary
