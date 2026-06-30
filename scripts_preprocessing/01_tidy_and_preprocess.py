"""
Full preprocessing pipeline for a single .iFLiP2 file.

Outputs to data/sessions/<session>/ (all files have no header — values only):

  data_params.csv            — header key/value pairs (with column labels)
  data_time.csv              — sample timestamps in seconds            (n_time x 1)
  data_intensity.csv         — total photon counts per timepoint       (n_time x 1)
  data_tcspc.csv             — TCSPC histogram matrix                  (n_time x n_lt_bins)
  data_marks.csv             — channel marks per timepoint             (n_time x 1)
  data_mpet.csv              — MPET in ns per timepoint                (n_time x 1)
  data_lt_singleExp_tau.csv  — per-timepoint tau1 from single-exp fit  (n_time x 1)
  data_lt_singleExp_t0.csv   — per-timepoint t0 from single-exp fit    (n_time x 1)
  data_lt_doubleExp_tau.csv  — per-timepoint tau1, tau2 from double-exp (n_time x 2)
  data_lt_doubleExp_t0.csv   — per-timepoint t0 from double-exp fit    (n_time x 1)

MPET and lifetime fits use background-corrected TCSPC data from the iFLiP2_Code
functions (calculate_mpet, fit_single_exp, fit_double_exp).

Each row of data_tcspc (one timepoint) is fitted independently. An accumulated
session-wide fit is run first and used as the warm-start beta0 for all
per-timepoint fits to accelerate convergence. Timepoints with fewer than
--min_counts photons are written as NaN.

Usage:
    python scripts_preprocessing/01_tidy_and_preprocess.py <input.iFLiP2> \\
        --spc_range 0.4 12.4 \\
        --t0 1.15 \\
        [--afterpulse_ratio 0.04] \\
        [--bg <background.iFLiP2>] \\
        [--pulse_interval 12.5] \\
        [--min_counts 10] \\
        [--sessions_dir data/sessions]
"""

import re
import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "functions"))
from read_iflip2 import read_iflip2
from calculate_mpet import calculate_mpet
from fit_lifetime import fit_single_exp, fit_double_exp

REPO_ROOT = Path(__file__).resolve().parent.parent


def session_name(stem):
    return re.sub(r"_\d+$", "", stem)


def flatten_header(d, prefix="header"):
    rows = []
    for k, v in d.items():
        key = f"{prefix}_{k}".lower()
        if isinstance(v, dict):
            if set(v.keys()) <= {"Value", "Text"}:
                leaf = v.get("Value", v.get("Text", ""))
                rows.append((key, leaf))
            else:
                rows.extend(flatten_header(v, prefix=key))
        else:
            rows.append((key, v))
    return rows


def run(args):
    filepath = Path(args.input)
    d = read_iflip2(filepath)
    data = d["data"]           # (n_lt_bins, n_time, n_channels), uint32
    header = d["header"]
    lt_time = d["LTTime"]      # ns, (n_lt_bins,)
    sample_time = d["sampleTime"]
    marks = d["marks"]
    n_time = data.shape[1]

    session = session_name(filepath.stem)
    out_dir = Path(args.sessions_dir) / session / "data_flipr"
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Background curve (optional) ---
    bg_lt_curve = None
    if args.bg:
        bg_d = read_iflip2(Path(args.bg))
        bg_lt_curve = bg_d["data"][:, :, 0].sum(axis=1).astype(float)

    # --- MPET (background-corrected, t0-subtracted) ---
    mpet, corrected_data = calculate_mpet(
        data, header,
        spc_range=(args.spc_range[0], args.spc_range[1]),
        t0=args.t0,
        bg_lt_curve=bg_lt_curve,
        after_pulse_ratio=args.afterpulse_ratio,
    )
    # mpet:           (n_time, n_channels)
    # corrected_data: (n_lt_bins, n_time, n_channels)

    # -----------------------------------------------------------------------
    # Session-wide fit on accumulated TCSPC — gives warm-start beta0 values
    # -----------------------------------------------------------------------
    lt_accumulated = corrected_data[:, :, 0].sum(axis=1)  # (n_lt_bins,)
    print(f"[{session}]  fitting session-wide TCSPC ...", flush=True)
    r1_session = fit_single_exp(lt_time, lt_accumulated, pulse_interval=args.pulse_interval)
    r2_session = fit_double_exp(lt_time, lt_accumulated, pulse_interval=args.pulse_interval)
    print(
        f"[{session}]  session fit: "
        f"single tau={r1_session['tau1']:.4f} ns | "
        f"double tau1={r2_session['tau1']:.4f} tau2={r2_session['tau2']:.4f} ns",
        flush=True,
    )

    # -----------------------------------------------------------------------
    # Per-timepoint lifetime fits
    # -----------------------------------------------------------------------
    tau_single  = np.full(n_time, np.nan)
    t0_single   = np.full(n_time, np.nan)
    tau1_double = np.full(n_time, np.nan)
    tau2_double = np.full(n_time, np.nan)
    t0_double   = np.full(n_time, np.nan)

    # Warm-start beta0 from session-wide fit (amplitudes will be re-scaled per timepoint)
    sig1 = r1_session["sigma"]
    sig2 = r2_session["sigma"]

    print(f"[{session}]  fitting {n_time} timepoints ...", flush=True)
    report_every = max(1, n_time // 10)

    for i in range(n_time):
        row = corrected_data[:, i, 0]   # (n_lt_bins,)
        total = float(row.sum())

        if total < args.min_counts:
            continue

        peak = max(float(row.max()), 1.0)

        # Single exp
        b0_1 = [peak, r1_session["tau1"], r1_session["t0"], sig1, 0.0]
        try:
            r = fit_single_exp(lt_time, row,
                               pulse_interval=args.pulse_interval, beta0=b0_1)
            tau_single[i] = r["tau1"]
            t0_single[i]  = r["t0"]
        except Exception:
            pass

        # Double exp
        b0_2 = [peak / 2, r2_session["tau1"],
                 peak / 2, r2_session["tau2"],
                 r2_session["t0"], sig2, 0.0]
        try:
            r = fit_double_exp(lt_time, row,
                               pulse_interval=args.pulse_interval, beta0=b0_2)
            tau1_double[i] = r["tau1"]
            tau2_double[i] = r["tau2"]
            t0_double[i]   = r["t0"]
        except Exception:
            pass

        if (i + 1) % report_every == 0:
            print(f"  {i + 1}/{n_time}", flush=True)

    # -----------------------------------------------------------------------
    # Write outputs — no header, values only
    # -----------------------------------------------------------------------

    # data_params.csv keeps column labels (key/value lookup table)
    param_rows = [("filename", filepath.stem)] + flatten_header(header)
    pd.DataFrame(param_rows, columns=["variable", "value"]).to_csv(
        out_dir / "data_params.csv", index=False
    )

    ch0 = data[:, :, 0]  # (n_lt_bins, n_timepoints)
    intensity = ch0.sum(axis=0).astype(float)

    pd.Series(sample_time).to_csv(
        out_dir / "data_time.csv", header=False, index=False)
    pd.Series(intensity).to_csv(
        out_dir / "data_intensity.csv", header=False, index=False)
    pd.Series(marks.astype(float)).to_csv(
        out_dir / "data_marks.csv", header=False, index=False)

    # TCSPC matrix: rows = timepoints, columns = LT bins
    pd.DataFrame(ch0.T).to_csv(
        out_dir / "data_tcspc.csv", header=False, index=False)

    # MPET: values only, one column per channel
    pd.DataFrame(mpet).to_csv(
        out_dir / "data_mpet.csv", header=False, index=False)

    # Per-timepoint lifetime fits
    pd.Series(tau_single).to_csv(
        out_dir / "data_lt_singleExp_tau.csv", header=False, index=False)
    pd.Series(t0_single).to_csv(
        out_dir / "data_lt_singleExp_t0.csv", header=False, index=False)
    pd.DataFrame({"tau1": tau1_double, "tau2": tau2_double}).to_csv(
        out_dir / "data_lt_doubleExp_tau.csv", header=False, index=False)
    pd.Series(t0_double).to_csv(
        out_dir / "data_lt_doubleExp_t0.csv", header=False, index=False)

    n_fit = int(np.isfinite(tau_single).sum())
    print(
        f"[{session}]  done — {n_fit}/{n_time} timepoints fitted | "
        f"MPET mean {np.nanmean(mpet):.4f} ns",
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Tidy and preprocess a .iFLiP2 file into session CSVs."
    )
    parser.add_argument("input", help=".iFLiP2 file path")
    parser.add_argument(
        "--spc_range", nargs=2, type=float, required=True,
        metavar=("T_MIN", "T_MAX"),
        help="SPC range in ns, e.g. --spc_range 0.4 12.4",
    )
    parser.add_argument(
        "--t0", type=float, required=True,
        help="Time offset to subtract from MPET, in ns",
    )
    parser.add_argument(
        "--afterpulse_ratio", type=float, default=0.0,
        help="After-pulse ratio (default: 0). Typically 0.04 (PDA44 setting 3) "
             "or 0.08 (setting 4)",
    )
    parser.add_argument(
        "--bg", default=None,
        help="Background .iFLiP2 file (all timepoints summed for background curve)",
    )
    parser.add_argument(
        "--pulse_interval", type=float, default=12.5,
        help="Laser repetition period in ns (default: 12.5 ns = 80 MHz)",
    )
    parser.add_argument(
        "--min_counts", type=float, default=10.0,
        help="Skip fitting for timepoints with fewer total photons (default: 10)",
    )
    parser.add_argument(
        "--sessions_dir",
        default=str(REPO_ROOT / "data" / "sessions"),
        help="Root sessions directory (default: data/sessions/)",
    )
    args = parser.parse_args()
    run(args)
