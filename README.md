# analysis-flipr

Analysis pipeline for **lifetime photometry (FLIPR)** data recorded with the iFLiP2 system.

---

## Repository structure

```
analysis-flipr/
├── functions/                  # Shared Python functions (ported from iFLiP2_Code)
│   ├── read_iflip2.py          # Parse .iFLiP2 binary files
│   ├── calculate_mpet.py       # MPET with dead-time + afterpulse correction
│   └── fit_lifetime.py         # Single- and double-exp lifetime fitting
├── scripts_preprocessing/      # Preprocessing pipeline scripts
│   └── 01_tidy_and_preprocess.py
├── scripts_analysis/           # Exploratory analysis, visualization, statistics
├── source/                     # Original iFLiP2_Code for Users (MATLAB reference)
├── data/                       # Working data — not tracked by git
│   └── sessions/
│       └── <session>/
│           └── data_flipr/     # All preprocessed outputs for this session
├── data_example/               # Small example datasets for testing
├── output/                     # Saved results: figures, tables, statistics
└── prompts/                    # AI prompt logs used during development
```

---

## Dependencies

```bash
pip install -r requirements.txt   # numpy, pandas, scipy
```

---

## Preprocessing

### `scripts_preprocessing/01_tidy_and_preprocess.py`

Full preprocessing pipeline for a single `.iFLiP2` session file. Reads the raw
binary, applies background and dead-time corrections, calculates MPET, and fits
per-timepoint lifetime decays. All outputs are written to
`data/sessions/<session>/data_flipr/` with **no header row**.

```bash
python scripts_preprocessing/01_tidy_and_preprocess.py <input.iFLiP2> \
    --spc_range 0.4 12.4 \
    --t0 1.15 \
    [--afterpulse_ratio 0.04] \
    [--bg <background.iFLiP2>] \
    [--pulse_interval 12.5] \
    [--min_counts 10] \
    [--sessions_dir data/sessions]
```

**Key parameters**

| Argument | Description |
|---|---|
| `--spc_range T_MIN T_MAX` | LT bins to include in MPET (ns). Typically `t0 − 0.6` to just before the edge artifact. |
| `--t0` | Time offset (ns) subtracted from MPET. Obtain from a session-wide single-exp fit. |
| `--afterpulse_ratio` | Detector afterpulse fraction (0.04 for PDA44 setting 3, 0.08 for setting 4). Default: 0. |
| `--bg` | Background `.iFLiP2` file; all timepoints are summed into a single background curve. |
| `--pulse_interval` | Laser repetition period in ns (default: 12.5 ns = 80 MHz). |
| `--min_counts` | Timepoints with fewer total photons are skipped and written as NaN (default: 10). |

**Outputs** (all in `data/sessions/<session>/data_flipr/`, no header)

| File | Shape | Description |
|---|---|---|
| `data_params.csv` | key/value | Flattened acquisition header (has column labels) |
| `data_time.csv` | n × 1 | Sample timestamps (s) |
| `data_intensity.csv` | n × 1 | Total photon counts per timepoint |
| `data_tcspc.csv` | n × 126 | Full TCSPC histogram matrix (one row per timepoint) |
| `data_marks.csv` | n × 1 | Channel event marks |
| `data_mpet.csv` | n × 1 | MPET per timepoint (ns, background-corrected, t0-subtracted) |
| `data_lt_singleExp_tau.csv` | n × 1 | Per-timepoint tau from single-exp fit (ns) |
| `data_lt_singleExp_t0.csv` | n × 1 | Per-timepoint t0 from single-exp fit (ns) |
| `data_lt_doubleExp_tau.csv` | n × 2 | Per-timepoint tau1, tau2 from double-exp fit (ns) |
| `data_lt_doubleExp_t0.csv` | n × 1 | Per-timepoint t0 from double-exp fit (ns) |

**Processing steps**

1. Parse `.iFLiP2` binary → raw TCSPC array `(n_lt_bins × n_time × n_channels)`
2. Background correction with dead-time and afterpulse models (`calculate_mpet`)
3. MPET: `sum(I(t)·t) / sum(I(t)) − t0` over the SPC range
4. Session-wide accumulated TCSPC fit → warm-start initial conditions
5. Per-timepoint single- and double-exponential fits (`fit_lifetime`)

All MPET and lifetime calculations use the ported iFLiP2_Code functions; no
ad-hoc implementations.

---

## Functions

Functions in `functions/` are ported directly from the MATLAB source in `source/iFLiP2_Code for Users/`.

| Module | Source | Description |
|---|---|---|
| `read_iflip2.py` | `iFLiP2_readData.m` | Parse `.iFLiP2` binary files |
| `calculate_mpet.py` | `iFLiP2_calculateMPET.m` | MPET with dead-time + afterpulse correction |
| `fit_lifetime.py` | `h_fitLifetimeBySingleExp.m` / `h_fitLifetimeByDoubleExp.m` | Exponential-Gaussian lifetime fitting |

---

## Contributing

Open issues or PRs on the [GitHub repository](https://github.com/gordon-laboratory/analysis-flipr).
