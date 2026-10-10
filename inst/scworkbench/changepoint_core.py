#!/usr/bin/env python3
"""Single-change-point numeric core for scwbtrajectory.

Phases (one process each, JSON file in / JSON file out):
  probe    - verify the ruptures API behaves as this operation uses it
  features - per-feature smoothing, one-breakpoint detection, delta, direction
  global   - cross-feature breakpoint histogram, density, threshold rank set
"""

import argparse
import hashlib
import json
import os
import sys

SCHEMA = "scwbtrajectory.changepoint_core.v1"

EXIT_IMPORT = 2
EXIT_API = 3
EXIT_REQUEST = 4
EXIT_COMPUTE = 5


class CoreError(Exception):
    def __init__(self, code, message):
        super().__init__(message)
        self.code = code
        self.message = message


def _core_sha256():
    digest = hashlib.sha256()
    with open(os.path.realpath(__file__), "rb") as handle:
        for block in iter(lambda: handle.read(1 << 16), b""):
            digest.update(block)
    return digest.hexdigest()


def runtime_identity():
    identity = {
        "argv0": sys.argv[0] if sys.argv else "",
        "executable": os.path.realpath(sys.executable),
        "executable_input": sys.executable,
        "python_version": "%d.%d.%d" % sys.version_info[:3],
        "changepoint_core_sha256": _core_sha256(),
    }
    try:
        from importlib import metadata
    except ImportError as error:
        raise CoreError(EXIT_IMPORT, "importlib.metadata unavailable: %s" % error)
    try:
        identity["ruptures_version"] = metadata.version("ruptures")
    except Exception:
        identity["ruptures_version"] = None
    try:
        import ruptures
        identity["ruptures_dunder_version"] = str(getattr(ruptures, "__version__", ""))
    except Exception:
        identity["ruptures_dunder_version"] = None
    for name in ("numpy", "scipy", "pandas"):
        try:
            identity["%s_version" % name] = metadata.version(name)
        except Exception:
            identity["%s_version" % name] = None
    return identity


def _window_class():
    try:
        from ruptures import Window
    except ImportError as error:
        raise CoreError(EXIT_IMPORT, "ruptures is not importable: %s" % error)
    except Exception as error:
        raise CoreError(EXIT_IMPORT, "ruptures import failed: %s" % error)
    return Window


def _window(width, jump, min_size):
    window_class = _window_class()
    try:
        return window_class(model="l1", width=int(width), jump=int(jump),
                            min_size=int(min_size))
    except Exception as error:
        raise CoreError(EXIT_API,
                        "Window(model='l1', width=%d, jump=%d, min_size=%d) "
                        "could not be constructed: %s" % (int(width), int(jump),
                                                          int(min_size), error))


def _predict_one(values, width, jump, min_size):
    import numpy as np
    window = _window(width, jump, min_size)
    try:
        raw = window.fit(np.asarray(values, dtype="float64")).predict(n_bkps=1)
    except Exception as error:
        raise CoreError(EXIT_API,
                        "predict(n_bkps=1) failed on a %d-point signal: %s"
                        % (len(values), error))
    try:
        bkps = [int(value) for value in raw]
    except Exception as error:
        raise CoreError(EXIT_API, "predict(n_bkps=1) returned non-integer %r" % (raw,))
    return bkps


def _to_vector(values, feature_id):
    """A JSON null is a missing observation, not a Python error.

    Null becomes NaN here so the caller's non_finite_policy decides what happens
    to the feature; a value that is neither null nor a number is a malformed
    request and is refused rather than silently turned into a missing value.
    """
    import numpy as np
    out = []
    missing = []
    for index, value in enumerate(values):
        if value is None:
            missing.append(index)
            out.append(float("nan"))
        elif isinstance(value, bool) or not isinstance(value, (int, float)):
            raise CoreError(EXIT_REQUEST,
                            "feature %r has a non-numeric value %r at position %d"
                            % (feature_id, value, index))
        else:
            out.append(float(value))
    return np.asarray(out, dtype="float64"), missing


def _window_params(params):
    width = int(params["window_width"])
    jump = max(1, int(params["window_jump"]))
    min_size = int(params["window_min_size"])
    if width < 2:
        raise CoreError(EXIT_REQUEST, "window_width must be >= 2")
    if min_size < 1:
        raise CoreError(EXIT_REQUEST, "window_min_size must be >= 1")
    if min_size >= width:
        raise CoreError(EXIT_REQUEST, "window_min_size must be < window_width")
    return width, jump, min_size


def phase_probe(request):
    params = request.get("params")
    if not isinstance(params, dict):
        raise CoreError(EXIT_REQUEST, "probe requires a params object")
    for name in ("window_width", "window_jump", "window_min_size"):
        if name not in params:
            raise CoreError(EXIT_REQUEST, "probe params is missing %r" % name)
    width, jump, min_size = _window_params(params)
    probe_n = max(4 * width, 200)
    planted = probe_n // 2
    import numpy as np
    signal = np.zeros(probe_n, dtype="float64")
    signal[planted:] = 1.0
    bkps = _predict_one(signal, width, jump, min_size)
    if len(bkps) != 2:
        raise CoreError(EXIT_API,
                        "probe expected the two-element [breakpoint, n_samples] "
                        "pair, got %r" % (bkps,))
    boundary, terminal = bkps[0], bkps[1]
    if terminal != probe_n:
        raise CoreError(EXIT_API,
                        "probe terminal sample count is %d, expected %d"
                        % (terminal, probe_n))
    if not 0 < boundary < probe_n:
        raise CoreError(EXIT_API,
                        "probe breakpoint %d is not interior to (0, %d)"
                        % (boundary, probe_n))
    if abs(boundary - planted) > jump:
        raise CoreError(EXIT_API,
                        "probe recovered boundary %d for the planted step at %d "
                        "(tolerance %d)" % (boundary, planted, jump))
    response = {
        "probe": {
            "probe_n": probe_n,
            "planted_boundary": planted,
            "recovered_boundary": boundary,
            "terminal_sample_count": terminal,
            "tolerance": jump,
            "signal": "deterministic_planted_step",
        }
    }
    response.update(runtime_identity())
    return response


def rolling_mean(values, window):
    window = int(window)
    if window <= 1:
        return [float(value) for value in values]
    try:
        import pandas as pd
    except ImportError as error:
        raise CoreError(EXIT_IMPORT, "pandas is not importable: %s" % error)
    series = pd.Series([float(value) for value in values], dtype="float64")
    # center=True with min_periods=0 reproduces the notebook's rmean /
    # get_smoothened behaviour: partial windows at both edges, never NaN.
    return series.rolling(window=window, center=True, min_periods=0).mean().tolist()


def _record(feature_id, b=None, delta=None, direction=None, analysed=False,
            reason=None, missing_values=0):
    return {"feature_id": feature_id, "b": b, "delta": delta,
            "direction": direction, "analysed": analysed, "reason": reason,
            "missing_values": int(missing_values)}


def phase_features(request):
    cells = request.get("cells")
    features = request.get("features")
    params = request.get("params")
    if not isinstance(cells, list) or not cells:
        raise CoreError(EXIT_REQUEST, "features phase requires a non-empty cells list")
    if not isinstance(features, list):
        raise CoreError(EXIT_REQUEST, "features phase requires a features list")
    if not isinstance(params, dict):
        raise CoreError(EXIT_REQUEST, "features phase requires a params object")
    for name in ("smooth_window", "window_width", "window_jump", "window_min_size",
                 "direction_window", "non_finite_policy", "degenerate_bkp_policy"):
        if name not in params:
            raise CoreError(EXIT_REQUEST, "features params is missing %r" % name)
    width, jump, min_size = _window_params(params)
    smooth_window = int(params["smooth_window"])
    direction_window = int(params["direction_window"])
    if direction_window < 1:
        raise CoreError(EXIT_REQUEST, "direction_window must be >= 1")
    non_finite_policy = str(params["non_finite_policy"])
    degenerate_policy = str(params["degenerate_bkp_policy"])
    if non_finite_policy not in ("error", "drop"):
        raise CoreError(EXIT_REQUEST, "non_finite_policy must be 'error' or 'drop'")
    if degenerate_policy not in ("error", "drop"):
        raise CoreError(EXIT_REQUEST, "degenerate_bkp_policy must be 'error' or 'drop'")
    import numpy as np
    count = len(cells)
    records = []
    for feature in features:
        if not isinstance(feature, dict) or "feature_id" not in feature:
            raise CoreError(EXIT_REQUEST, "every feature needs a feature_id")
        feature_id = str(feature["feature_id"])
        values = feature.get("values")
        if not isinstance(values, list) or len(values) != count:
            raise CoreError(EXIT_REQUEST,
                            "feature %r has %s values for %d cells"
                            % (feature_id, len(values) if isinstance(values, list)
                               else "non-list", count))
        vector, missing = _to_vector(values, feature_id)
        if not bool(np.all(np.isfinite(vector))):
            if non_finite_policy == "error":
                other = [int(index) for index in np.flatnonzero(~np.isfinite(vector))
                         if index not in set(missing)]
                positions = sorted(missing + other)
                shown = ", ".join(str(index) for index in positions[:8])
                if len(positions) > 8:
                    shown += ", ... (%d more)" % (len(positions) - 8)
                raise CoreError(EXIT_COMPUTE,
                                "feature %r has %d missing (null) and %d non-finite "
                                "(NaN/Inf) values out of %d cells, at positions [%s]; "
                                "non_finite_policy='error' refuses the run "
                                "(non_finite_policy='drop' publishes the feature as "
                                "unanalysed)" % (feature_id, len(missing), len(other),
                                                 count, shown))
            records.append(_record(feature_id, reason="non_finite",
                                   missing_values=len(missing)))
            continue
        smoothed = rolling_mean(vector, smooth_window)
        if not all(value == value for value in smoothed):
            raise CoreError(EXIT_COMPUTE,
                            "smoothing produced NaN for feature %r" % feature_id)
        bkps = _predict_one(smoothed, width, jump, min_size)
        boundary = bkps[0] if bkps else count
        if not bkps or boundary >= count:
            if degenerate_policy == "error":
                raise CoreError(EXIT_COMPUTE,
                                "feature %r produced no interior breakpoint "
                                "(degenerate_bkp_policy='error'); predict returned %r"
                                % (feature_id, bkps))
            records.append(_record(feature_id, reason="degenerate"))
            continue
        array = np.asarray(smoothed, dtype="float64")
        after = array[boundary:boundary + direction_window]
        before = array[max(0, boundary - direction_window):boundary]
        if after.size == 0 or before.size == 0:
            if degenerate_policy == "error":
                raise CoreError(EXIT_COMPUTE,
                                "feature %r has an empty direction window at boundary %d"
                                % (feature_id, boundary))
            records.append(_record(feature_id, reason="degenerate"))
            continue
        delta = float(after.mean() - before.mean())
        records.append(_record(feature_id, b=int(boundary), delta=delta,
                               direction="up" if delta > 0 else "down",
                               analysed=True))
    return {"records": records}


def _density(counts, count_length, smooth_window, scale):
    """Notebook convention: slot j of the returned vector carries boundary j+1.

    counts buckets are 0..count_length (length count_length + 1); buckets 0 and
    count_length are dropped, the remainder is smoothed then scaled, and a
    literal 0 is appended for the final cell slot.
    """
    if len(counts) != count_length + 1:
        raise CoreError(EXIT_REQUEST,
                        "count array has %d entries, expected %d"
                        % (len(counts), count_length + 1))
    inner = [float(value) for value in counts[1:count_length]]
    raw = [scale * value for value in rolling_mean(inner, smooth_window)] if inner else []
    if len(raw) != max(0, count_length - 1):
        raise CoreError(EXIT_COMPUTE,
                        "raw density has %d entries, expected %d"
                        % (len(raw), max(0, count_length - 1)))
    return raw + [0.0]


def phase_global(request):
    cell_ids = request.get("cell_ids")
    params = request.get("params")
    counts = request.get("counts")
    if not isinstance(cell_ids, list) or not cell_ids:
        raise CoreError(EXIT_REQUEST, "global phase requires a non-empty cell_ids list")
    if not isinstance(params, dict):
        raise CoreError(EXIT_REQUEST, "global phase requires a params object")
    for name in ("density_smooth_window", "density_scale", "change_threshold"):
        if name not in params:
            raise CoreError(EXIT_REQUEST, "global params is missing %r" % name)
    if not isinstance(counts, dict):
        raise CoreError(EXIT_REQUEST, "global phase requires a counts object")
    n = len(cell_ids)
    smooth_window = int(params["density_smooth_window"])
    scale = float(params["density_scale"])
    threshold = float(params["change_threshold"])
    raw = {}
    for key in ("all", "up", "down"):
        if key not in counts:
            raise CoreError(EXIT_REQUEST, "counts is missing %r" % key)
        values = counts[key]
        if not isinstance(values, list):
            raise CoreError(EXIT_REQUEST, "counts[%r] must be a list" % key)
        raw[key] = _density(values, n, smooth_window, scale)
    # D1: publish boundary b at cell rank b instead of the notebook's b-1.
    aligned = {}
    for key in ("all", "up", "down"):
        shifted = [0.0] + raw[key][0:n - 1]
        if len(shifted) != n:
            raise CoreError(EXIT_COMPUTE,
                            "aligned density for %r has %d entries, expected %d"
                            % (key, len(shifted), n))
        aligned[key] = shifted
    # A boundary rank is significant for a direction when that direction's own
    # density reaches the threshold there. The combined ("all") density mixes the
    # two directions, so a strong up population would otherwise make a lone down
    # breakpoint look significant; it stays published as a diagnostic only.
    significant = {}
    for key in ("all", "up", "down"):
        significant[key] = [index for index in range(n)
                            if aligned[key][index] >= threshold]

    def peak(values):
        best = 0
        for index in range(1, len(values)):
            if values[index] > values[best]:
                best = index
        return best

    return {
        "density_all": aligned["all"],
        "density_up": aligned["up"],
        "density_down": aligned["down"],
        "significant_ranks": significant["all"],
        "significant_ranks_up": significant["up"],
        "significant_ranks_down": significant["down"],
        "significance": {
            "region_selection": "direction_matched_density_v1",
            "rule": "rank r is significant for direction d when the aligned "
                    "density of d at r is >= change_threshold",
            "combined_density_role": "diagnostic_only",
            "threshold_units": "absolute units of a density_scale-scaled "
                               "breakpoint count; not scale free",
        },
        "alignment": {
            "published": "aligned_boundary_at_cell_rank",
            "notebook_raw": "slot_j_carries_boundary_j_plus_1",
            "shift": "aligned[j] = raw[j - 1]; aligned[0] = 0",
            "raw_peak_rank_all": peak(raw["all"]),
            "aligned_peak_rank_all": peak(aligned["all"]),
            "raw_peak_rank_up": peak(raw["up"]),
            "aligned_peak_rank_up": peak(aligned["up"]),
            "raw_peak_rank_down": peak(raw["down"]),
            "aligned_peak_rank_down": peak(aligned["down"]),
        },
    }


PHASES = {"probe": phase_probe, "features": phase_features, "global": phase_global}


def run(phase, request):
    response = {"schema": SCHEMA, "phase": phase, "ok": True}
    response.update(PHASES[phase](request))
    return response


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", required=True, choices=sorted(PHASES))
    parser.add_argument("--request", required=True)
    parser.add_argument("--response", required=True)
    args = parser.parse_args(argv)
    try:
        if args.phase not in PHASES:
            raise CoreError(EXIT_REQUEST, "unknown phase %r" % args.phase)
        with open(args.request, "r", encoding="utf-8") as handle:
            request = json.load(handle)
        if not isinstance(request, dict):
            raise CoreError(EXIT_REQUEST, "request is not an object")
        if request.get("schema") != SCHEMA:
            raise CoreError(EXIT_REQUEST,
                            "request schema %r is not %r" % (request.get("schema"), SCHEMA))
        if request.get("phase") != args.phase:
            raise CoreError(EXIT_REQUEST,
                            "request phase %r disagrees with --phase %r"
                            % (request.get("phase"), args.phase))
        response = run(args.phase, request)
    except CoreError as error:
        response = {"schema": SCHEMA, "phase": args.phase, "ok": False,
                    "error": error.code, "message": error.message}
    except Exception as error:
        response = {"schema": SCHEMA, "phase": args.phase, "ok": False,
                    "error": EXIT_COMPUTE, "message": "%s: %s" % (type(error).__name__, error)}
    with open(args.response, "w", encoding="utf-8") as handle:
        json.dump(response, handle)
    return 0 if response.get("ok") else response.get("error", EXIT_COMPUTE)


if __name__ == "__main__":
    sys.exit(main())
