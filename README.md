# analysis-flipr

Analysis pipeline for **lifetime photometry (FLIPR)** data.

## Repository structure

```
analysis-flipr/
├── scripts_preprocessing/   # Preprocessing pipelines (signal extraction, normalization, etc.)
├── scripts_analysis/        # Exploratory analysis, visualization, and statistics
├── data_example/            # Example datasets for testing and demonstration
├── data/                    # Working data (not tracked by git)
├── output/                  # Saved results: figures, tables, statistics
└── prompts/                 # AI prompt logs used during development
```

## Getting started

1. Clone the repo:
   ```bash
   git clone https://github.com/gordon-laboratory/analysis-flipr.git
   ```
2. Install dependencies (environment file TBD):
   ```bash
   pip install -r requirements.txt
   ```
3. Run preprocessing on your data:
   ```bash
   python scripts_preprocessing/<script>.py --input data/<your_file>
   ```

## Contributing
Open issues or PRs on the [GitHub repository](https://github.com/gordon-laboratory/analysis-flipr).
