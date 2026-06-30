import re
import numpy as np
from pathlib import Path


def read_iflip2(filepath):
    """
    Read a .iFLiP2 lifetime photometry file.

    Returns a dict with keys:
        header      - nested dict of all header fields
        data        - np.ndarray, shape (LTCurveLength, timepoints, nChannels), uint32
        LTTime      - 1D array of lifetime-axis time points in ns
        sampleTime  - 1D array of sample timestamps in seconds
        marks       - 1D array of channel marks (one value per timepoint)
    """
    filepath = Path(filepath)
    with open(filepath, "rb") as fh:
        header_lines = []
        while True:
            line = fh.readline().decode("ascii", errors="replace").rstrip("\r\n")
            assert line != "", "Reached end of file before 'header_end' marker"
            if line == "header_end":
                break
            header_lines.append(line)
        raw_data = np.frombuffer(fh.read(), dtype=np.uint32)

    header = _parse_header(header_lines)

    n_channels = int(header["acq"]["nChannels"])
    lt_curve_len = int(header["acq"]["LTCurveLength"]) + 1  # +1 for marks row
    lt_resolution = float(header["acq"]["LTResolution"])
    sampling_freq = float(header["state"]["samplingFreq"]["Value"])

    n_timepoints = raw_data.size / n_channels / lt_curve_len
    if n_timepoints != round(n_timepoints):
        raise ValueError(f"Data size mismatch: {raw_data.size} values cannot be "
                         f"divided evenly into [{lt_curve_len} x {n_channels} x N]")
    n_timepoints = int(n_timepoints)

    # MATLAB reshape is column-major; replicate with order='F' then transpose axes
    # MATLAB: permute(reshape(data, [LTCurveLength, nChannels, dataTimeLength]), [1 3 2])
    # → Python axes after permute: [LTCurve, time, channel]
    data = (raw_data
            .reshape((lt_curve_len, n_channels, n_timepoints), order="F")
            .transpose(0, 2, 1))  # [LTCurveLength, timepoints, nChannels]

    lt_time = np.arange(lt_curve_len - 1) * lt_resolution          # ns
    sample_time = (np.arange(n_timepoints) + 0.5) / sampling_freq  # seconds

    marks = data[-1, :, 0].copy()   # last row, all timepoints, ch0 (same across channels)
    data = data[:-1, :, :]          # strip marks row

    return {
        "header": header,
        "data": data,
        "LTTime": lt_time,
        "sampleTime": sample_time,
        "marks": marks,
    }


# ---------------------------------------------------------------------------
# Header parser
# ---------------------------------------------------------------------------

def _parse_header(lines):
    """Parse MATLAB-style struct assignment lines into a nested dict."""
    header = {}
    for line in lines:
        line = line.strip()
        if not line:
            continue
        # Match: header.a.b.c... = <value>  OR  header.a.b.c....Text = <value>
        m = re.match(r"^header\.(.+?)\s*=\s*(.+);$", line)
        if not m:
            continue
        key_path = m.group(1).split(".")
        raw_val = m.group(2).strip()
        value = _parse_value(raw_val)
        _set_nested(header, key_path, value)
    return header


def _parse_value(raw):
    """Convert a MATLAB literal string to a Python scalar or string."""
    if raw.startswith("'") and raw.endswith("'"):
        return raw[1:-1]
    try:
        return int(raw)
    except ValueError:
        pass
    try:
        return float(raw)
    except ValueError:
        pass
    return raw  # fall back to raw string


def _set_nested(d, keys, value):
    for key in keys[:-1]:
        d = d.setdefault(key, {})
    d[keys[-1]] = value
