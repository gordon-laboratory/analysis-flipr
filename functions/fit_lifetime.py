"""
Port of h_fitLifetimeBySingleExp.m and h_fitLifetimeByDoubleExp.m.

Fits a TCSPC lifetime decay curve (accumulated over many time-points) to a
single or double exponential model convolved with a Gaussian IRF, with
optional pre-pulse contributions summed from previous laser periods.

Model (single exp):
    y(x) = sum_{k=0}^{n} amp1/2 * exp(sigma^2/(2*tau1^2) - (x - t0 + k*T)/tau1)
                                 * erfc((sigma^2 - tau1*(x - t0 + k*T)) / (sqrt(2)*tau1*sigma))
           + bg

where T = pulse_interval (ns) and n = ceil(10 * tau1 / T).

The double-exp model adds a second exponential component with (amp2, tau2).

Reference:  h_fitLifetimeBySingleExp.m / h_fitLifetimeByDoubleExp.m
            (iFLiP2 software, Bhattacharya / Bhattacharyya labs)
"""

import numpy as np
from scipy.optimize import least_squares
from scipy.special import erfc


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def fit_single_exp(lt_time, intensity, pulse_interval=12.5,
                   beta0=None, is_fix=None):
    """
    Fit a single-exponential + Gaussian-IRF model to a TCSPC decay curve.

    Parameters
    ----------
    lt_time : array_like, shape (n,)
        TCSPC time axis in ns (e.g., from read_iflip2 LTTime).
    intensity : array_like, shape (n,)
        Photon counts at each LT bin.
    pulse_interval : float, optional
        Laser repetition period in ns (1 / sync_rate * 1e9). Default: 12.5 ns
        (80 MHz).
    beta0 : array_like, shape (5,), optional
        Initial parameters [amp1, tau1, t0, sigma, bg]. Auto-estimated if None.
    is_fix : array_like of bool/int, shape (5,), optional
        Which parameters to hold fixed at their beta0 values.
        Order: [amp1, tau1, t0, sigma, bg].
        Default: [0, 0, 0, 0, 1] — fix bg at 0.

    Returns
    -------
    result : dict
        amp1, tau1, t0, sigma, bg, converge (bool), fitted_x, fitted_y,
        original_x, original_y, method.
    """
    lt_time = np.asarray(lt_time, dtype=float)
    intensity = np.asarray(intensity, dtype=float)

    if beta0 is None:
        beta0 = _init_single(lt_time, intensity)
    if is_fix is None:
        is_fix = [0, 0, 0, 0, 1]  # fix bg at 0

    beta0 = np.asarray(beta0, dtype=float)
    is_fix = np.asarray(is_fix, dtype=bool)

    model_fn = lambda beta, x: _model_single(x, *beta, pulse_interval)
    betahat, converge = _nlinfit(lt_time, intensity, model_fn, beta0, is_fix)

    # Evaluate on 10× oversampled grid for plotting
    t_step = float(np.mean(np.diff(lt_time)))
    x_fine = np.arange(lt_time[0], lt_time[-1], 0.1 * t_step)
    y_fine = _model_single(x_fine, *betahat, pulse_interval)

    return {
        "amp1": betahat[0], "tau1": betahat[1],
        "t0": betahat[2], "sigma": betahat[3], "bg": betahat[4],
        "beta": betahat, "converge": converge,
        "fitted_x": x_fine, "fitted_y": y_fine,
        "original_x": lt_time, "original_y": intensity,
        "method": "Single Exp",
    }


def fit_double_exp(lt_time, intensity, pulse_interval=12.5,
                   beta0=None, is_fix=None):
    """
    Fit a double-exponential + Gaussian-IRF model to a TCSPC decay curve.

    Parameters
    ----------
    lt_time : array_like, shape (n,)
        TCSPC time axis in ns.
    intensity : array_like, shape (n,)
        Photon counts at each LT bin.
    pulse_interval : float, optional
        Laser repetition period in ns. Default: 12.5 ns (80 MHz).
    beta0 : array_like, shape (7,), optional
        Initial parameters [amp1, tau1, amp2, tau2, t0, sigma, bg].
        Auto-estimated if None.
    is_fix : array_like of bool/int, shape (7,), optional
        Which parameters to hold fixed at their beta0 values.
        Order: [amp1, tau1, amp2, tau2, t0, sigma, bg].
        Default: [0, 0, 0, 0, 0, 0, 1] — fix bg at 0.

    Returns
    -------
    result : dict
        amp1, tau1, amp2, tau2, t0, sigma, bg, converge (bool),
        fitted_x, fitted_y, original_x, original_y, method.
    """
    lt_time = np.asarray(lt_time, dtype=float)
    intensity = np.asarray(intensity, dtype=float)

    if beta0 is None:
        beta0 = _init_double(lt_time, intensity)
    if is_fix is None:
        is_fix = [0, 0, 0, 0, 0, 0, 1]  # fix bg at 0

    beta0 = np.asarray(beta0, dtype=float)
    is_fix = np.asarray(is_fix, dtype=bool)

    model_fn = lambda beta, x: _model_double(x, *beta, pulse_interval)
    betahat, converge = _nlinfit(lt_time, intensity, model_fn, beta0, is_fix)

    t_step = float(np.mean(np.diff(lt_time)))
    x_fine = np.arange(lt_time[0], lt_time[-1], 0.1 * t_step)
    y_fine = _model_double(x_fine, *betahat, pulse_interval)

    return {
        "amp1": betahat[0], "tau1": betahat[1],
        "amp2": betahat[2], "tau2": betahat[3],
        "t0": betahat[4], "sigma": betahat[5], "bg": betahat[6],
        "beta": betahat, "converge": converge,
        "fitted_x": x_fine, "fitted_y": y_fine,
        "original_x": lt_time, "original_y": intensity,
        "method": "Double Exp",
    }


# ---------------------------------------------------------------------------
# Model functions
# ---------------------------------------------------------------------------

def _exp_gauss_component(x, amp, tau, t0, sigma, pulse_interval):
    """
    Single exponential convolved with Gaussian IRF, summed over pre-pulses.

    Returns contribution for one (amp, tau) component.
    x must be a 1-D array.
    """
    n = max(2, round(tau / pulse_interval * 10))
    k = np.arange(n + 1)
    # t_offsets[k] = -t0 + k * T  →  (x + t_offset) = x - t0 + k*T
    t_offsets = -t0 + k * pulse_interval          # (n+1,)
    u = x[:, np.newaxis] + t_offsets[np.newaxis, :]  # (n_pts, n+1)

    y1 = amp * np.exp(sigma**2 / (2 * tau**2) - u / tau)
    y2 = erfc((sigma**2 - tau * u) / (np.sqrt(2) * tau * sigma))
    return (y1 * y2).sum(axis=1) / 2


def _model_single(x, amp1, tau1, t0, sigma, bg, pulse_interval):
    return _exp_gauss_component(x, amp1, tau1, t0, sigma, pulse_interval) + bg


def _model_double(x, amp1, tau1, amp2, tau2, t0, sigma, bg, pulse_interval):
    c1 = _exp_gauss_component(x, amp1, tau1, t0, sigma, pulse_interval)
    c2 = _exp_gauss_component(x, amp2, tau2, t0, sigma, pulse_interval)
    return c1 + c2 + bg


# ---------------------------------------------------------------------------
# Auto-initialisation helpers
# ---------------------------------------------------------------------------

def _init_single(t, y):
    """Auto-estimate beta0 = [amp1, tau1, t0, sigma, bg] for single exp."""
    peak_idx = int(np.argmax(y))
    estimated_t0 = t[peak_idx] - 0.5
    mpet = np.sum(t * y) / np.sum(y) - estimated_t0 if np.sum(y) > 0 else 1.0
    return [float(np.max(y)), max(mpet, 0.1), estimated_t0, 0.13, 0.0]


def _init_double(t, y):
    """Auto-estimate beta0 = [amp1, tau1, amp2, tau2, t0, sigma, bg] for double exp."""
    peak_idx = int(np.argmax(y))
    max_y = float(np.max(y))
    estimated_t0 = t[peak_idx] - 0.5
    mpet = np.sum(t * y) / np.sum(y) - estimated_t0 if np.sum(y) > 0 else 1.0
    mpet = max(mpet, 0.1)
    return [max_y / 2, mpet * 2, max_y / 2, mpet / 2, estimated_t0, 0.13, 0.0]


# ---------------------------------------------------------------------------
# Optimizer
# ---------------------------------------------------------------------------

def _nlinfit(x_data, y_data, model_fn, beta0, is_fix):
    """
    Weighted nonlinear least-squares fit via scipy.optimize.least_squares.

    Replicates the weighting from ss_nlinfit:
        weight = max(sqrt(intensity), 1) / sqrt(max(intensity))
    Residuals = (y - model) / weight  (Poisson chi-weighting).

    Fixed parameters (is_fix=True) are held at their beta0 values;
    only free parameters are passed to the optimizer.

    Returns
    -------
    betahat : ndarray, shape (n_params,)
        Full parameter vector (fixed params restored).
    converge : bool
    """
    beta0 = np.asarray(beta0, dtype=float)
    is_fix = np.asarray(is_fix, dtype=bool)
    free_mask = ~is_fix

    # Poisson weighting: higher intensity → higher weight → smaller relative residual
    # Clip to zero before sqrt — background-subtracted rows can be slightly negative.
    y_clipped = np.clip(y_data, 0.0, None)
    max_y = max(float(y_clipped.max()), 1.0)
    weight = np.maximum(np.sqrt(y_clipped), 1.0) / np.sqrt(max_y)

    def full_model(free_params):
        full = beta0.copy()
        full[free_mask] = free_params
        return model_fn(full, x_data)

    def residual(free_params):
        return (y_data - full_model(free_params)) / weight

    x0 = beta0[free_mask]
    opt = least_squares(residual, x0, method="lm", max_nfev=10_000)

    betahat = beta0.copy()
    betahat[free_mask] = opt.x
    converge = opt.success

    return betahat, converge
