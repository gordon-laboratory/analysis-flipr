"""
Port of iFLiP2_calculateMPET.m.

Calculate Mean Photon Emission Time (MPET) from raw TCSPC lifetime histograms,
with dead-time and afterpulse background corrections.

Formula:  MPET = sum(I(t) * t) / sum(I(t)) - t0
"""

import numpy as np


def calculate_mpet(data, header, spc_range, t0,
                   bg_lt_curve=None, after_pulse_ratio=0.0):
    """
    Calculate MPET from raw lifetime curves.

    Parameters
    ----------
    data : ndarray, shape (n_lt_bins, n_time, n_channels), uint32
        Raw TCSPC histograms as returned by read_iflip2.
    header : dict
        Header dict as returned by read_iflip2.
    spc_range : (float, float)
        (t_min, t_max) in ns — range of LT bins to include. Typically
        t_min ~ t0 - 0.6, t_max as large as possible before the edge artifact.
    t0 : float
        Time offset (ns) to subtract from MPET. Obtain from software fit or
        from h_fitLifetimeBySingleExp.
    bg_lt_curve : ndarray, shape (n_lt_bins,), optional
        Background TCSPC histogram (single time-point). Default: zeros.
    after_pulse_ratio : float, optional
        After-pulse ratio (typically 0.04 for PDA44 setting 3, 0.08 for
        setting 4). Default: 0 (no after-pulse correction).

    Returns
    -------
    mpet : ndarray, shape (n_time, n_channels), float64
        Mean photon emission time minus t0, in ns.
        NaN where total photon count in the SPC range is 0.
    corrected_data : ndarray, shape (n_lt_bins, n_time, n_channels), float64
        Background-subtracted TCSPC histograms.
    """
    n_lt_bins = data.shape[0]

    sampling_freq = float(header["state"]["samplingFreq"]["Value"])
    init = header.get("init", {})
    dead_time = float(init["deadTime_ns"]) / 1e9 if "deadTime_ns" in init else 25e-9

    lt_resolution = float(header["acq"]["LTResolution"])
    lt_time = np.arange(n_lt_bins) * lt_resolution  # ns

    if bg_lt_curve is None:
        bg_lt_curve = np.zeros(n_lt_bins, dtype=np.float64)

    corrected_data = _calc_effective_bg(
        data, bg_lt_curve, dead_time, sampling_freq, after_pulse_ratio
    )

    idx = (lt_time >= spc_range[0]) & (lt_time <= spc_range[1])
    if not idx.any():
        raise ValueError(
            f"No LTTime bins fall within spc_range {spc_range}. "
            f"LTTime spans [{lt_time[0]:.3f}, {lt_time[-1]:.3f}] ns."
        )

    lt_time_use = lt_time[idx]                          # (n_sel,)
    data_use = corrected_data[idx, :, :]                # (n_sel, n_time, n_channels)

    total_counts = data_use.sum(axis=0)                 # (n_time, n_channels)
    weighted_sum = (
        data_use * lt_time_use[:, np.newaxis, np.newaxis]
    ).sum(axis=0)                                       # (n_time, n_channels)

    mpet = np.full(total_counts.shape, np.nan, dtype=np.float64)
    valid = total_counts > 0
    mpet[valid] = weighted_sum[valid] / total_counts[valid] - t0

    return mpet, corrected_data


# ---------------------------------------------------------------------------
# Internal background-correction helper
# ---------------------------------------------------------------------------

def _calc_effective_bg(lt_histogram, measured_bg, dead_time, sampling_freq,
                       afterpulse_ratio):
    """
    Apply dead-time and afterpulse corrections to background.

    Implements internal_calcEffectiveBG from iFLiP2_calculateMPET.m.

    Dead-time model:  count_efficiency = 1 - raw_rate * dead_time
    Afterpulse model: flat contribution proportional to signal intensity.

    Parameters
    ----------
    lt_histogram : ndarray, shape (n_lt_bins, n_time, n_channels), uint32
    measured_bg  : ndarray, shape (n_lt_bins,)
        Measured background TCSPC histogram (single time-point, no corrections).
    dead_time    : float, seconds
    sampling_freq : float, Hz
    afterpulse_ratio : float

    Returns
    -------
    corrected_data : ndarray, same shape as lt_histogram, float64
    """
    lt_histogram = lt_histogram.astype(np.float64)
    n_lt_bins = lt_histogram.shape[0]

    # --- Background dead-time correction ---
    # E = 1 - r0 * dt  (first-order approximation)
    bg_count_efficiency = 1.0 - measured_bg.sum() * sampling_freq * dead_time
    true_bg = measured_bg / bg_count_efficiency  # (n_lt_bins,)

    # --- Per-timepoint dead-time correction ---
    raw_intensity_per_sec = lt_histogram.sum(axis=0) * sampling_freq  # (n_time, n_channels)
    count_efficiency = 1.0 - raw_intensity_per_sec * dead_time        # (n_time, n_channels)

    # Effective background: scale true_bg by per-timepoint count efficiency
    # Broadcasting: (n_lt_bins, 1, 1) * (1, n_time, n_channels)
    effective_bg = (
        true_bg[:, np.newaxis, np.newaxis]
        * count_efficiency[np.newaxis, :, :]
    )  # (n_lt_bins, n_time, n_channels)

    # Effective BG integrated intensity (counts per second)
    effective_bg_int = true_bg.sum() * count_efficiency * sampling_freq  # (n_time, n_channels)

    # Afterpulse: flat across LT bins, proportional to net signal intensity
    afterpulse_bg = (
        (raw_intensity_per_sec - effective_bg_int) * afterpulse_ratio
        / (sampling_freq * n_lt_bins)
    )  # (n_time, n_channels)

    total_bg = effective_bg + afterpulse_bg[np.newaxis, :, :]
    return lt_histogram - total_bg
