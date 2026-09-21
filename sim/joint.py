#!/usr/bin/env python3
"""Joint hybrid–smooth toy: occult mode guards on a Lyapunov collocation.

The smooth proliferative field is a planar tumour–effector equation.
The collocation follows the meshless orbital-derivative construction of
Argáez, Giesl and Hafstein (Wendland kernel, speed-normalised field,
a failing set where the derivative will not stay negative). That failing
set is a collocation defect on a fixed window. It is not a computer-assisted
proof and not a certified Conley set.

Quiescence, angiogenic pause, and immune-held latency are discrete modes.
Their guards are functions of the shared constants. No mode adds a
coordinate, and no mode introduces a new numerical parameter.

The run is deterministic. The label 20260921 is stored and is not an
input to the arithmetic.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve
from scipy.ndimage import label
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
RUN_LABEL = 20260921

# Shared constants. Time is scaled by the tumour growth rate.
# T is burden in units of carrying capacity. E is effector density
# in units of a reference density. Not estimated from a patient or a cell line.
A_PAR = 1.0
B_PAR = 0.5
C_PAR = 0.65
D_PAR = 0.2

T_STAR = D_PAR * B_PAR / (C_PAR * A_PAR - D_PAR)
E_STAR = (1.0 - T_STAR) * (B_PAR + T_STAR) / A_PAR
Y_DET = 0.15

T_LO, T_HI = 0.03, 0.58
E_LO, E_HI = 0.25, 0.92
N_PRIMARY = 34
N_COARSE = 22
N_FINE = 42
SUPPORT_MULT = 9.0
GAMMA_MENU = (-0.55, -0.40, -0.25)
STENCIL_M = 4
STENCIL_STEP = 0.4
DELTA2 = 1e-8
N_BINARY_ITERS = 5
GUARD_SAMPLES = 401

# Classification radii, in units of the larger grid spacing.
INTERIOR_NEAR = 1.0
INTERIOR_FAR = 2.0
INTERFACE_NEAR = 1.5
TRANSIENT_FAR = 2.0


def field_p(z: np.ndarray) -> np.ndarray:
    """Proliferative mode. The only mode with a positive coexistence state."""
    z = np.asarray(z, dtype=float)
    t = z[..., 0]
    e = z[..., 1]
    kill = A_PAR * t * e / (B_PAR + t)
    dT = t * (1.0 - t) - kill
    dE = C_PAR * kill - D_PAR * e
    return np.stack([dT, dE], axis=-1)


def field_q(z: np.ndarray) -> np.ndarray:
    """Quiescence: growth removed, recruitment removed. Kill and clearance remain."""
    z = np.asarray(z, dtype=float)
    t = z[..., 0]
    e = z[..., 1]
    kill = A_PAR * t * e / (B_PAR + t)
    return np.stack([-kill, -D_PAR * e], axis=-1)


def field_a(z: np.ndarray) -> np.ndarray:
    """Angiogenic pause: tumour coordinate held, effector clearance unchanged."""
    z = np.asarray(z, dtype=float)
    e = z[..., 1]
    dT = np.zeros_like(e)
    return np.stack([dT, -D_PAR * e], axis=-1)


def field_i(z: np.ndarray) -> np.ndarray:
    """Immune-held structure: growth removed, recruitment kept, clearance kept."""
    z = np.asarray(z, dtype=float)
    t = z[..., 0]
    e = z[..., 1]
    kill = A_PAR * t * e / (B_PAR + t)
    return np.stack([-kill, C_PAR * kill - D_PAR * e], axis=-1)


FIELDS = {"P": field_p, "Q": field_q, "A": field_a, "I": field_i}
GUARD_PARTNER = {"cyc": "Q", "imm": "I", "ang": "A"}


def jacobian_p(z: np.ndarray) -> np.ndarray:
    t, e = float(z[0]), float(z[1])
    den = B_PAR + t
    dfT = (1.0 - 2.0 * t) - A_PAR * e * B_PAR / den**2
    dfE = -A_PAR * t / den
    dgT = C_PAR * A_PAR * e * B_PAR / den**2
    dgE = C_PAR * A_PAR * t / den - D_PAR
    return np.array([[dfT, dfE], [dgT, dgE]], dtype=float)


def phi_derivatives(s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """First and second derivatives of unscaled Wendland psi_{4,2}."""
    s = np.asarray(s, dtype=float)
    d1 = np.zeros_like(s)
    d2 = np.zeros_like(s)
    mask = (s > 0.0) & (s < 1.0)
    u = 1.0 - s[mask]
    d1[mask] = -(u**5 / 5.0 - 11.0 * u**6 / 30.0 + u**7 / 6.0)
    d2[mask] = u**4 - (11.0 / 5.0) * u**5 + (7.0 / 6.0) * u**6
    return d1, d2


def psi_pair(dist: np.ndarray, shape: float) -> tuple[np.ndarray, np.ndarray]:
    dist = np.asarray(dist, dtype=float)
    p1 = np.zeros_like(dist)
    p2 = np.zeros_like(dist)
    zero = dist <= 1e-14
    p1[zero] = -(shape**2) / 30.0
    p2[zero] = shape**4
    positive = ~zero
    radius = np.zeros_like(dist)
    radius[positive] = shape * dist[positive]
    inside = positive & (radius < 1.0)
    if np.any(inside):
        sm = radius[inside]
        d1, d2 = phi_derivatives(sm)
        p1[inside] = shape**2 * d1 / sm
        p2[inside] = shape**4 * (sm * d2 - d1) / sm**3
    return p1, p2


def normalise(f: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(f, axis=-1)
    return f / np.sqrt(DELTA2 + norm**2)[..., None]


def assemble(points: np.ndarray, f_hat: np.ndarray, shape: float) -> np.ndarray:
    diff = points[:, None, :] - points[None, :, :]
    dist = np.linalg.norm(diff, axis=-1)
    p1, p2 = psi_pair(dist, shape)
    fi = f_hat[:, None, :]
    fj = f_hat[None, :, :]
    dot_i = np.sum(diff * fi, axis=-1)
    dot_j = np.sum((-diff) * fj, axis=-1)
    dot_f = np.sum(fi * fj, axis=-1)
    return p2 * dot_i * dot_j - p1 * dot_f


def values_and_derivative(
    query: np.ndarray,
    f_query: np.ndarray,
    points: np.ndarray,
    f_hat: np.ndarray,
    beta: np.ndarray,
    shape: float,
) -> tuple[np.ndarray, np.ndarray]:
    diff = query[:, None, :] - points[None, :, :]
    dist = np.linalg.norm(diff, axis=-1)
    p1, p2 = psi_pair(dist, shape)
    fj = f_hat[None, :, :]
    fq = f_query[:, None, :]
    orbital = -p1 * np.sum(fq * fj, axis=-1) + p2 * np.sum(diff * fq, axis=-1) * np.sum(
        (-diff) * fj, axis=-1
    )
    coef = np.sum((-diff) * fj, axis=-1) * p1
    return coef @ beta, orbital @ beta


def cartesian_grid(n: int) -> tuple[np.ndarray, float, np.ndarray, np.ndarray]:
    ts = np.linspace(T_LO, T_HI, n)
    es = np.linspace(E_LO, E_HI, n)
    grid = np.stack(np.meshgrid(ts, es, indexing="xy"), axis=-1).reshape(-1, 2)
    spacing = float(max(ts[1] - ts[0], es[1] - es[0]))
    return grid, spacing, ts, es


def stencil_blocks(points: np.ndarray, f_hat: np.ndarray, spacing: float) -> np.ndarray:
    speed = np.linalg.norm(f_hat, axis=1)
    direction = np.zeros_like(f_hat)
    moving = speed > 1e-14
    direction[moving] = f_hat[moving] / speed[moving, None]
    step = STENCIL_STEP * spacing
    blocks = []
    for k in range(1, STENCIL_M + 1):
        blocks.append(points + k * step * direction)
        blocks.append(points - k * step * direction)
    return np.stack(blocks, axis=1)


def stencil_means(
    blocks: np.ndarray,
    points: np.ndarray,
    f_hat: np.ndarray,
    beta: np.ndarray,
    shape: float,
) -> np.ndarray:
    sample = blocks.reshape(-1, 2)
    _, orbital = values_and_derivative(sample, normalise(field_p(sample)), points, f_hat, beta, shape)
    return orbital.reshape(blocks.shape[0], blocks.shape[1]).mean(axis=1)


def sample_orbit(eq: np.ndarray) -> tuple[np.ndarray, float, dict]:
    def ode(_t, z):
        return field_p(np.asarray(z, dtype=float))

    sol = solve_ivp(
        ode,
        (0.0, 420.0),
        eq * np.array([1.12, 0.88]),
        rtol=1e-8,
        atol=1e-8,
        dense_output=True,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    times = np.linspace(300.0, 420.0, 9000)
    path = sol.sol(times)
    series = path[0]
    crossings = np.where((series[:-1] < eq[0]) & (series[1:] >= eq[0]))[0]
    if len(crossings) < 3:
        raise RuntimeError("periodic orbit was not detected by the section")
    periods = np.diff(times[crossings])
    period = float(np.median(periods[-4:]))
    i0, i1 = crossings[-3], crossings[-2]
    orbit = path[:, i0 : i1 + 1].T
    meta = {
        "n_samples": int(len(orbit)),
        "T_min": float(np.min(orbit[:, 0])),
        "T_max": float(np.max(orbit[:, 0])),
        "E_min": float(np.min(orbit[:, 1])),
        "E_max": float(np.max(orbit[:, 1])),
        "period_samples_std": float(np.std(periods[-4:])),
        "section": "upward crossings of T = T*",
    }
    return orbit, period, meta


def self_test() -> None:
    rng = np.random.default_rng(0)
    centres = rng.normal(size=(10, 2)) * 0.2
    beta = rng.normal(size=10)
    shape = 1.7

    def f_test(x):
        return np.stack([-x[..., 1], x[..., 0]], axis=-1)

    f_centres = f_test(centres)
    x = np.array([[0.15, -0.08]])
    f = f_test(x)
    _, orbital = values_and_derivative(x, f, centres, f_centres, beta, shape)
    eps = 1e-6
    forward, _ = values_and_derivative(x + eps * f, f_test(x + eps * f), centres, f_centres, beta, shape)
    backward, _ = values_and_derivative(x - eps * f, f_test(x - eps * f), centres, f_centres, beta, shape)
    numeric = (forward - backward) / (2.0 * eps)
    rel = abs(float(orbital[0] - numeric[0])) / max(1.0, abs(float(numeric[0])))
    if rel > 1e-6:
        raise RuntimeError(f"orbital derivative disagrees with a finite difference ({rel})")
    pts = np.array([[0.2, 0.3], [0.25, 0.45], [0.4, 0.35], [0.33, 0.55], [0.5, 0.4]])
    fh = normalise(field_p(pts))
    matrix = assemble(pts, fh, 2.0)
    if np.max(np.abs(matrix - matrix.T)) > 1e-12:
        raise RuntimeError("collocation matrix is not symmetric")
    coeff = solve(matrix, -np.ones(len(pts)), assume_a="pos")
    _, got = values_and_derivative(pts, fh, pts, fh, coeff, 2.0)
    if np.max(np.abs(got + 1.0)) > 1e-8:
        raise RuntimeError("collocation residual at the nodes is not zero")


def orbit_flatness(beta, points, f_hat, shape, orbit) -> dict:
    sample = orbit[:: max(1, len(orbit) // 160)]
    values, orbital = values_and_derivative(sample, normalise(field_p(sample)), points, f_hat, beta, shape)
    return {
        "V_peak_to_peak": float(np.max(values) - np.min(values)),
        "V_net": float(values[-1] - values[0]),
        "V_mean": float(np.mean(values)),
        "orbital_median": float(np.median(orbital)),
        "orbital_p10": float(np.quantile(orbital, 0.1)),
        "orbital_p90": float(np.quantile(orbital, 0.9)),
        "orbital_max": float(np.max(orbital)),
        "fraction_positive": float(np.mean(orbital > 0.0)),
        "sample_T": sample[:, 0],
        "sample_V": values,
        "sample_orbital": orbital,
    }


def connected_components(fail: np.ndarray, n: int, points: np.ndarray, orbit: np.ndarray, eq: np.ndarray) -> list[dict]:
    mask = fail.reshape(n, n)
    labeled, nlab = label(mask, structure=np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]]))
    flat = labeled.reshape(-1)
    orbit_tree = cKDTree(orbit[::2])
    comps = []
    for k in range(1, nlab + 1):
        idx = np.where(flat == k)[0]
        pts = points[idx]
        d_orb = orbit_tree.query(pts)[0]
        d_eq = np.linalg.norm(pts - eq, axis=1)
        med_o = float(np.median(d_orb))
        med_f = float(np.median(d_eq))
        if med_o + 1e-15 < med_f:
            kind = "orbit"
        elif med_f + 1e-15 < med_o:
            kind = "focus"
        else:
            kind = "ambiguous"
        comps.append(
            {
                "id": int(k),
                "size": int(len(idx)),
                "kind": kind,
                "median_dist_orbit": med_o,
                "median_dist_focus": med_f,
                "centroid": [float(np.mean(pts[:, 0])), float(np.mean(pts[:, 1]))],
                "_idx": idx,
            }
        )
    comps.sort(key=lambda c: -c["size"])
    return comps


def classify_points(query: np.ndarray, comps: list[dict], points: np.ndarray, spacing: float) -> list[dict]:
    query = np.atleast_2d(np.asarray(query, dtype=float))
    trees = []
    for comp in comps:
        trees.append(cKDTree(points[comp["_idx"]]))
    rows = []
    for p in query:
        dists = [float(tree.query(p)[0]) for tree in trees] if trees else []
        if not dists:
            rows.append(
                {
                    "T": float(p[0]),
                    "E": float(p[1]),
                    "class": "transient",
                    "component_id": None,
                    "component_kind": None,
                    "dist_nearest": None,
                    "dist_second": None,
                }
            )
            continue
        order = np.argsort(dists)
        d1 = dists[order[0]]
        d2 = dists[order[1]] if len(order) > 1 else float("inf")
        comp = comps[int(order[0])]
        comp2 = comps[int(order[1])] if len(order) > 1 else None
        if d1 <= INTERFACE_NEAR * spacing and d2 <= INTERFACE_NEAR * spacing and comp["kind"] != (comp2 or {}).get("kind"):
            klass = "interface"
            cid, kind = None, None
        elif d1 <= INTERIOR_NEAR * spacing and d2 > INTERIOR_FAR * spacing:
            klass = "interior"
            cid, kind = comp["id"], comp["kind"]
        elif d1 > TRANSIENT_FAR * spacing:
            klass = "transient"
            cid, kind = None, None
        else:
            klass = "margin"
            cid, kind = comp["id"], comp["kind"]
        rows.append(
            {
                "T": float(p[0]),
                "E": float(p[1]),
                "class": klass,
                "component_id": cid,
                "component_kind": kind,
                "dist_nearest": d1,
                "dist_second": None if not np.isfinite(d2) else d2,
            }
        )
    return rows


def summarise_classes(rows: list[dict], orbit: np.ndarray | None = None) -> dict:
    n = len(rows)
    keys = ("interior", "interface", "transient", "margin")
    counts = {k: int(sum(r["class"] == k for r in rows)) for k in keys}
    interior_kinds = {}
    interior_ids = {}
    for r in rows:
        if r["class"] == "interior":
            interior_kinds[r["component_kind"]] = interior_kinds.get(r["component_kind"], 0) + 1
            interior_ids[str(r["component_id"])] = interior_ids.get(str(r["component_id"]), 0) + 1
    out = {
        "n": n,
        "counts": counts,
        "fractions": {k: (counts[k] / n if n else None) for k in keys},
        "interior_by_kind": interior_kinds,
        "interior_by_component_id": interior_ids,
    }
    if orbit is not None and n:
        tree = cKDTree(orbit[::2])
        pts = np.array([[r["T"], r["E"]] for r in rows])
        dist = tree.query(pts)[0]
        interior = np.array([r["class"] == "interior" for r in rows])
        out["dist_to_orbit_median"] = float(np.median(dist))
        out["dist_to_orbit_min"] = float(np.min(dist))
        if np.any(interior):
            out["interior_dist_to_orbit_min"] = float(np.min(dist[interior]))
            out["interior_dist_to_orbit_median"] = float(np.median(dist[interior]))
            out["interior_dist_to_orbit_max"] = float(np.max(dist[interior]))
        else:
            out["interior_dist_to_orbit_min"] = None
            out["interior_dist_to_orbit_median"] = None
            out["interior_dist_to_orbit_max"] = None
    return out


def guard_segment(name: str) -> np.ndarray:
    if name == "cyc":
        e = np.linspace(E_LO, E_HI, GUARD_SAMPLES)
        return np.column_stack([np.full(GUARD_SAMPLES, B_PAR), e])
    if name == "ang":
        e = np.linspace(E_LO, E_HI, GUARD_SAMPLES)
        return np.column_stack([np.full(GUARD_SAMPLES, T_STAR), e])
    if name == "imm":
        t = np.linspace(T_LO, T_HI, GUARD_SAMPLES)
        return np.column_stack([t, np.full(GUARD_SAMPLES, E_STAR)])
    raise KeyError(name)


def observable(name: str, z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    if name == "cyc":
        return z[..., 0] - B_PAR
    if name == "ang":
        return z[..., 0] - T_STAR
    if name == "imm":
        return z[..., 1] - E_STAR
    raise KeyError(name)


def orbit_crossings(orbit: np.ndarray, name: str) -> np.ndarray:
    h = observable(name, orbit)
    hits = []
    for i in range(len(orbit) - 1):
        if h[i] == 0.0:
            hits.append(orbit[i])
            continue
        if h[i] * h[i + 1] < 0.0:
            w = abs(h[i]) / (abs(h[i]) + abs(h[i + 1]))
            hits.append((1.0 - w) * orbit[i] + w * orbit[i + 1])
    if len(hits) == 0:
        return np.zeros((0, 2))
    arr = np.vstack(hits)
    # Drop a duplicate closing sample if the stored segment repeats its first point.
    if len(arr) >= 2 and np.linalg.norm(arr[0] - arr[-1]) < 1e-6:
        arr = arr[:-1]
    # Unique within a small tolerance, preserving order.
    kept = [arr[0]]
    for p in arr[1:]:
        if np.linalg.norm(p - kept[-1]) > 1e-4:
            kept.append(p)
    return np.vstack(kept)


def filippov_report(name: str, partner: str) -> dict:
    pts = guard_segment(name)
    if name in ("cyc", "ang"):
        normal = np.array([1.0, 0.0])
    else:
        normal = np.array([0.0, 1.0])
    fp = field_p(pts)
    fq = FIELDS[partner](pts)
    sp = fp @ normal
    sq = fq @ normal
    opposite = sp * sq < 0.0
    same = sp * sq > 0.0
    either_zero = (np.abs(sp) <= 1e-10) | (np.abs(sq) <= 1e-10)
    # Sliding vector when a convex weight kills the normal component.
    # lam * sp + (1-lam) * sq = 0 => lam = sq / (sq - sp) when signs oppose.
    lam = np.full(len(pts), np.nan)
    slide_speed = np.full(len(pts), np.nan)
    tangent = np.array([-normal[1], normal[0]])
    good = opposite & ~either_zero
    lam[good] = sq[good] / (sq[good] - sp[good])
    mix = lam[good, None] * fp[good] + (1.0 - lam[good, None]) * fq[good]
    slide_speed[good] = mix @ tangent
    return {
        "partner_mode": partner,
        "normal": normal.tolist(),
        "n": int(len(pts)),
        "fraction_opposite": float(np.mean(opposite)),
        "fraction_same_sign": float(np.mean(same)),
        "fraction_either_normal_near_zero": float(np.mean(either_zero)),
        "normal_P_median_abs": float(np.median(np.abs(sp))),
        "normal_partner_median_abs": float(np.median(np.abs(sq))),
        "sliding_speed_median_abs": float(np.nanmedian(np.abs(slide_speed))) if np.any(good) else None,
        "sliding_speed_max_abs": float(np.nanmax(np.abs(slide_speed))) if np.any(good) else None,
        "fraction_sliding_speed_below_1e-3": float(np.mean(np.abs(slide_speed[good]) < 1e-3)) if np.any(good) else None,
    }


def attractiveness(eq: np.ndarray, orbit: np.ndarray, period: float) -> dict:
    """Distance to the sampled orbit after several periods, from inside and outside."""

    def ode(_t, z):
        return field_p(np.asarray(z, dtype=float))

    tree = cKDTree(orbit[::2])
    rows = []
    starts = {
        "inside_offset": eq + np.array([0.04, -0.03]),
        "outside_corner": np.array([0.55, 0.80]),
    }
    for name, z0 in starts.items():
        sol = solve_ivp(ode, (0.0, 8.0 * period), z0, rtol=1e-7, atol=1e-7, dense_output=True)
        times = np.linspace(0.0, 8.0 * period, 9)
        path = sol.sol(times).T
        dists = tree.query(path)[0]
        rows.append(
            {
                "name": name,
                "start": z0.tolist(),
                "dist_start": float(dists[0]),
                "dist_end": float(dists[-1]),
                "decreased": bool(dists[-1] < dists[0]),
            }
        )
    return {"trials": rows}


def hybrid_switches(orbit: np.ndarray, period: float, eq: np.ndarray) -> list[dict]:
    """Cross a declared guard in mode P, then follow the partner field."""

    def events_for(active: str):
        specs = []
        for name in ("cyc", "imm", "ang"):
            def ev(_t, z, name=name):
                return float(observable(name, np.asarray(z)))

            ev.terminal = True
            ev.direction = 0.0
            specs.append((name, ev))
        return specs

    out = []
    # One start on the orbit, nudged off the angiogenic section so t=0 is not a hit.
    # The stored orbit begins at an upward crossing of T = T*.
    i_mid = len(orbit) // 5
    starts = [
        ("orbit_sample", orbit[i_mid], "P"),
        ("outside_corner", np.array([0.55, 0.80]), "P"),
        ("near_focus", eq + np.array([0.05, 0.04]), "P"),
    ]
    for label, z0, mode in starts:
        names_events = events_for(mode)

        def ode(_t, z, mode=mode):
            return FIELDS[mode](np.asarray(z, dtype=float))

        sol = solve_ivp(
            ode,
            (0.0, 3.0 * period),
            np.asarray(z0, dtype=float),
            events=[ev for _, ev in names_events],
            rtol=1e-7,
            atol=1e-7,
            dense_output=False,
        )
        hit_name = None
        hit_state = None
        hit_time = None
        for (name, _ev), te, ye in zip(names_events, sol.t_events, sol.y_events):
            if len(te) == 0:
                continue
            if hit_time is None or float(te[0]) < hit_time:
                # Ignore a hit that is numerically the initial point.
                if float(te[0]) < 1e-6:
                    continue
                hit_time = float(te[0])
                hit_name = name
                hit_state = np.asarray(ye[0], dtype=float)
        row = {
            "label": label,
            "start": np.asarray(z0, dtype=float).tolist(),
            "start_mode": mode,
            "hit": hit_name,
            "hit_time": hit_time,
            "hit_state": None if hit_state is None else hit_state.tolist(),
        }
        if hit_state is not None:
            partner = GUARD_PARTNER[hit_name]

            def ode_partner(_t, z, partner=partner):
                return FIELDS[partner](np.asarray(z, dtype=float))

            sol2 = solve_ivp(
                ode_partner,
                (0.0, 12.0),
                hit_state,
                rtol=1e-7,
                atol=1e-7,
                dense_output=True,
            )
            ts = np.linspace(0.0, 12.0, 241)
            path = sol2.sol(ts).T
            t_series = path[:, 0]
            below = np.where(t_series < Y_DET)[0]
            row["partner_mode"] = partner
            row["state_after_12"] = path[-1].tolist()
            row["distance_from_hit_after_12"] = float(np.linalg.norm(path[-1] - hit_state))
            row["min_T_after_switch"] = float(np.min(t_series))
            row["occult_within_12"] = bool(len(below) > 0)
            row["time_to_detection_floor"] = None if len(below) == 0 else float(ts[below[0]])
            row["_path_pause"] = path
            row["_hit"] = hit_state
        out.append(row)
    return out


def rk4_batch(z: np.ndarray, dt: float) -> np.ndarray:
    k1 = field_p(z)
    k2 = field_p(z + 0.5 * dt * k1)
    k3 = field_p(z + 0.5 * dt * k2)
    k4 = field_p(z + dt * k3)
    return z + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def return_diagnostic(points: np.ndarray, fail: np.ndarray, period: float) -> dict:
    dt = 0.05
    nsteps = int(np.ceil(2.0 * period / dt))
    z = points.copy()
    max_exc = np.zeros(len(z))
    close = np.zeros(len(z), dtype=bool)
    left = np.zeros(len(z), dtype=bool)
    t = 0.0
    lo_t, hi_t = T_LO - 0.02, T_HI + 0.02
    lo_e, hi_e = E_LO - 0.02, E_HI + 0.02
    for _ in range(nsteps):
        z = rk4_batch(z, dt)
        t += dt
        dist = np.linalg.norm(z - points, axis=1)
        max_exc = np.maximum(max_exc, dist)
        if 0.85 * period <= t <= 1.15 * period:
            close |= dist < 0.06
        left |= (z[:, 0] < lo_t) | (z[:, 0] > hi_t) | (z[:, 1] < lo_e) | (z[:, 1] > hi_e)
    returned = (max_exc >= 0.12) & close & ~left
    def rate(mask):
        if not np.any(mask):
            return None
        return float(np.mean(returned[mask]))
    return {
        "dt": dt,
        "steps": nsteps,
        "excursion_floor": 0.12,
        "return_radius": 0.06,
        "return_window": [0.85, 1.15],
        "fail_return_fraction": rate(fail),
        "pass_return_fraction": rate(~fail),
        "fail_left_fraction": float(np.mean(left[fail])) if np.any(fail) else None,
        "pass_left_fraction": float(np.mean(left[~fail])) if np.any(~fail) else None,
        "note": "A finite-time return is not an epsilon-chain. Leaving the window is not chain recurrence.",
    }


def collocate(n: int, gamma: float, orbit: np.ndarray, eq: np.ndarray, beta0=None, pack=None) -> dict:
    if pack is None:
        points, spacing, ts, es = cartesian_grid(n)
        shape = 1.0 / (SUPPORT_MULT * spacing)
        f_hat = normalise(field_p(points))
        matrix = assemble(points, f_hat, shape)
        symmetry = float(np.max(np.abs(matrix - matrix.T)))
        eig_min = float(np.min(np.linalg.eigvalsh(matrix)))
        cond = float(np.linalg.cond(matrix))
        beta0 = solve(matrix, -np.ones(len(points)), assume_a="pos")
        blocks = stencil_blocks(points, f_hat, spacing)
        score0 = stencil_means(blocks, points, f_hat, beta0, shape)
        pack = {
            "points": points,
            "spacing": spacing,
            "ts": ts,
            "es": es,
            "shape": shape,
            "f_hat": f_hat,
            "matrix": matrix,
            "symmetry": symmetry,
            "eig_min": eig_min,
            "cond": cond,
            "beta0": beta0,
            "blocks": blocks,
            "score0": score0,
        }
    fail = pack["score0"] > gamma
    rhs = np.where(fail, 0.0, -1.0)
    beta = solve(pack["matrix"], rhs, assume_a="pos")
    node_orbital = pack["matrix"] @ beta
    flat = orbit_flatness(beta, pack["points"], pack["f_hat"], pack["shape"], orbit)
    comps = connected_components(fail, n, pack["points"], orbit, eq)
    cover = None
    if np.any(fail):
        cover = float(np.mean(cKDTree(pack["points"][fail]).query(orbit[::2])[0] <= 2.0 * pack["spacing"]))
    else:
        cover = 0.0
    public_comps = []
    for c in comps:
        public_comps.append({k: v for k, v in c.items() if k != "_idx"})
    return {
        "n": n,
        "n_points": int(len(pack["points"])),
        "spacing": pack["spacing"],
        "shape": pack["shape"],
        "support_radius": float(1.0 / pack["shape"]),
        "matrix_symmetry_max": pack["symmetry"],
        "condition_number": pack["cond"],
        "min_eigenvalue": pack["eig_min"],
        "gamma": gamma,
        "fail_count": int(np.sum(fail)),
        "fail_fraction": float(np.mean(fail)),
        "orbit_cover_within_2_spacing": cover,
        "node_residual_inf": float(np.max(np.abs(node_orbital - rhs))),
        "flatness": {k: v for k, v in flat.items() if not k.startswith("sample_")},
        "components": public_comps,
        "_pack": pack,
        "_fail": fail,
        "_beta": beta,
        "_flat_samples": flat,
        "_comps": comps,
    }


def select_gamma(records: list[dict]) -> dict:
    eligible = [r for r in records if r["orbit_cover_within_2_spacing"] >= 0.95]
    pool = eligible if eligible else records
    pool = sorted(pool, key=lambda r: (r["flatness"]["V_peak_to_peak"], -r["orbit_cover_within_2_spacing"]))
    return pool[0]


def transient_drop(beta, points, f_hat, shape, period: float) -> dict:
    def ode(_t, z):
        return field_p(np.asarray(z, dtype=float))

    start = np.array([0.55, 0.80])
    horizon = 8.0 * period
    sol = solve_ivp(ode, (0.0, horizon), start, rtol=1e-7, atol=1e-7, dense_output=True)
    times = np.linspace(0.0, horizon, 800)
    path = sol.sol(times).T
    values = np.empty(len(path))
    for i in range(0, len(path), 100):
        sl = slice(i, i + 100)
        values[sl], _ = values_and_derivative(
            path[sl], normalise(field_p(path[sl])), points, f_hat, beta, shape
        )
    i_per = int(np.searchsorted(times, period))
    return {
        "start": start.tolist(),
        "V_start": float(values[0]),
        "V_after_one_period": float(values[i_per]),
        "drop_one_period": float(values[0] - values[i_per]),
        "V_end": float(values[-1]),
        "drop_total": float(values[0] - values[-1]),
        "times": times,
        "values": values,
        "path": path,
    }


def axial_equilibria() -> dict:
    rows = {}
    for name, z in {
        "origin": np.array([0.0, 0.0]),
        "carrying": np.array([1.0, 0.0]),
        "coexistence": np.array([T_STAR, E_STAR]),
    }.items():
        jac = jacobian_p(z)
        ev = np.linalg.eigvals(jac)
        rows[name] = {
            "state": z.tolist(),
            "field_norm": float(np.linalg.norm(field_p(z))),
            "trace": float(np.trace(jac)),
            "det": float(np.linalg.det(jac)),
            "eigenvalues_real": [float(np.real(v)) for v in ev],
            "eigenvalues_imag": [float(np.imag(v)) for v in ev],
        }
    return rows


def pause_equilibria_in_window() -> dict:
    """Scan a fine mesh for near-zeros of the three pause fields. Structural, not a proof."""
    ts = np.linspace(T_LO, T_HI, 80)
    es = np.linspace(E_LO, E_HI, 80)
    grid = np.stack(np.meshgrid(ts, es, indexing="xy"), axis=-1).reshape(-1, 2)
    out = {}
    for mode, fn in (("Q", field_q), ("A", field_a), ("I", field_i)):
        speed = np.linalg.norm(fn(grid), axis=1)
        i = int(np.argmin(speed))
        out[mode] = {
            "min_speed_on_mesh": float(speed[i]),
            "argmin_state": grid[i].tolist(),
            "positive_equilibrium_expected": False,
        }
    return out


def jsonable(obj):
    if isinstance(obj, dict):
        return {k: jsonable(v) for k, v in obj.items() if not str(k).startswith("_")}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def draw(primary, orbit, eq, crossings, hybrid, drop) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    points = primary["_pack"]["points"]
    fail = primary["_fail"]
    comps = primary["_comps"]
    spacing = primary["spacing"]
    n = primary["n"]
    kind_of = np.full(len(points), "", dtype=object)
    for comp in comps:
        kind_of[comp["_idx"]] = comp["kind"]

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 11,
        }
    )

    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.scatter(points[~fail, 0], points[~fail, 1], s=6, c="#d9d9d9", linewidths=0, label="passing nodes", zorder=1)
    styles = {"orbit": ("#1f4e79", "orbit defect"), "focus": ("#8c2f39", "focus defect"), "ambiguous": ("#6b4c9a", "ambiguous defect")}
    for kind, (color, label) in styles.items():
        m = fail & (kind_of == kind)
        if np.any(m):
            ax.scatter(points[m, 0], points[m, 1], s=12, c=color, linewidths=0, label=label, zorder=2)
    ax.plot(orbit[:, 0], orbit[:, 1], color="black", lw=1.2, label="periodic orbit", zorder=3)
    ax.scatter([eq[0]], [eq[1]], c="black", marker="*", s=80, label="coexistence focus", zorder=4)
    guard_style = {
        "cyc": ("#2a9d8f", "cycling guard T = b"),
        "imm": ("#e09f3e", "immune guard E = E*"),
        "ang": ("#6a4c93", "angiogenic guard T = T*"),
    }
    for name, (color, label) in guard_style.items():
        seg = guard_segment(name)
        ax.plot(seg[:, 0], seg[:, 1], color=color, lw=1.3, ls="--", label=label, zorder=3)
        if name in crossings and len(crossings[name]):
            ax.scatter(
                crossings[name][:, 0],
                crossings[name][:, 1],
                s=36,
                facecolors="none",
                edgecolors=color,
                linewidths=1.4,
                zorder=5,
            )
    ax.set_xlim(T_LO, T_HI)
    ax.set_ylim(E_LO, E_HI)
    ax.set_xlabel("tumour burden T")
    ax.set_ylabel("effector density E")
    ax.set_title("Collocation defect with hybrid guards")
    ax.legend(loc="upper right", fontsize=7.5, frameon=True)
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    fig.savefig(FIG / "partition_overlay.png", dpi=150)
    plt.close(fig)

    # Class stacks.
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    names = ["cyc", "imm", "ang"]
    classes = ["interior", "interface", "margin", "transient"]
    colors = {"interior": "#1f4e79", "interface": "#e09f3e", "margin": "#7a7a7a", "transient": "#c5c5c5"}
    x = np.arange(len(names))
    bottom = np.zeros(len(names))
    summaries = primary["_guard_summary"]
    for klass in classes:
        vals = np.array([summaries[name]["fractions"][klass] for name in names])
        ax.bar(x, vals, bottom=bottom, color=colors[klass], width=0.62, label=klass)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(["cycling T = b", "immune E = E*", "angiogenic T = T*"])
    ax.set_ylabel("fraction of guard samples")
    ax.set_ylim(0, 1.02)
    ax.set_title("Where each guard sits on the defect partition")
    ax.legend(fontsize=8, frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.18))
    fig.tight_layout()
    fig.savefig(FIG / "guard_classes.png", dpi=150)
    plt.close(fig)

    flat = primary["_flat_samples"]
    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    ax.axhline(0.0, color="#888", lw=0.8)
    ax.plot(flat["sample_T"], flat["sample_orbital"], color="#1f4e79", lw=1.1)
    ax.set_xlabel("tumour burden along the sampled orbit")
    ax.set_ylabel("rebuilt orbital derivative")
    ax.set_title("Orbital derivative along the orbit after one rebuild")
    fig.tight_layout()
    fig.savefig(FIG / "orbital_derivative.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    ax.plot(flat["sample_T"], flat["sample_V"], color="black", lw=1.2, label="on the orbit")
    # transient against time, separate axis note: plot V against arc progress
    ax.set_xlabel("tumour burden along the sampled orbit")
    ax.set_ylabel("rebuilt V")
    ax.set_title("Rebuilt function along the periodic orbit")
    fig.tight_layout()
    fig.savefig(FIG / "orbit_level.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    ax.plot(drop["times"] / primary["_period"], drop["values"], color="#8c2f39", lw=1.2)
    ax.set_xlabel("time / period")
    ax.set_ylabel("rebuilt V")
    ax.set_title("Transient from (0.55, 0.80), mode P")
    fig.tight_layout()
    fig.savefig(FIG / "transient_decrease.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.scatter(points[~fail, 0], points[~fail, 1], s=6, c="#eeeeee", linewidths=0, zorder=1)
    for kind, (color, label) in styles.items():
        m = fail & (kind_of == kind)
        if np.any(m):
            ax.scatter(points[m, 0], points[m, 1], s=10, c=color, linewidths=0, label=label, zorder=2)
    ax.plot(orbit[:, 0], orbit[:, 1], color="black", lw=1.0, label="periodic orbit", zorder=3)
    colors_h = {"orbit_sample": "#e09f3e", "outside_corner": "#2a9d8f", "near_focus": "#9b2226"}
    for row in hybrid:
        if row.get("_hit") is None:
            continue
        path = row["_path_pause"]
        color = colors_h.get(row["label"], "black")
        ax.plot(path[:, 0], path[:, 1], color=color, lw=1.3, label=f"after {row['label']} → {row['partner_mode']}", zorder=4)
        ax.scatter([row["_hit"][0]], [row["_hit"][1]], c=color, s=40, zorder=5)
    ax.set_xlim(T_LO - 0.01, T_HI + 0.01)
    ax.set_ylim(E_LO - 0.02, E_HI + 0.02)
    ax.set_xlabel("tumour burden T")
    ax.set_ylabel("effector density E")
    ax.set_title("Pause fields after a guard hit")
    ax.legend(loc="best", fontsize=7.5)
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    fig.savefig(FIG / "hybrid_traces.png", dpi=150)
    plt.close(fig)

    # Silence unused in case of lint.
    del n, spacing


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    self_test()
    eq = np.array([T_STAR, E_STAR], dtype=float)
    if abs(float(np.linalg.norm(field_p(eq)))) > 1e-12:
        raise RuntimeError("closed-form coexistence state is not an equilibrium")
    orbit, period, orbit_meta = sample_orbit(eq)
    attract = attractiveness(eq, orbit, period)
    equilibria = axial_equilibria()
    pause_eq = pause_equilibria_in_window()

    by_grid = {}
    chosen = {}
    for n in (N_COARSE, N_PRIMARY, N_FINE):
        pack = None
        records = []
        first = None
        for gamma in GAMMA_MENU:
            rec = collocate(n, gamma, orbit, eq, pack=pack)
            pack = rec["_pack"]
            if first is None:
                first = rec
            records.append(rec)
        winner = select_gamma(records)
        by_grid[n] = records
        chosen[n] = winner

    primary = chosen[N_PRIMARY]
    # Binary rebuild drift on the selected cut. Not the reported function.
    pack = primary["_pack"]
    beta_it = pack["beta0"].copy()
    history = []
    fail_it = primary["_fail"]
    for it in range(N_BINARY_ITERS):
        sc = stencil_means(pack["blocks"], pack["points"], pack["f_hat"], beta_it, pack["shape"])
        fail_it = sc > primary["gamma"]
        beta_it = solve(pack["matrix"], np.where(fail_it, 0.0, -1.0), assume_a="pos")
        comps_it = connected_components(fail_it, N_PRIMARY, pack["points"], orbit, eq)
        history.append(
            {
                "iteration": it + 1,
                "fail_fraction": float(np.mean(fail_it)),
                "n_components": len(comps_it),
                "largest": int(comps_it[0]["size"]) if comps_it else 0,
            }
        )

    crossings = {}
    crossing_rows = {}
    guard_rows = {}
    guard_summary = {}
    filippov = {}
    for name in ("cyc", "imm", "ang"):
        pts = guard_segment(name)
        rows = classify_points(pts, primary["_comps"], pack["points"], primary["spacing"])
        guard_rows[name] = rows
        guard_summary[name] = summarise_classes(rows, orbit)
        hits = orbit_crossings(orbit, name)
        crossings[name] = hits
        hit_class = classify_points(hits, primary["_comps"], pack["points"], primary["spacing"]) if len(hits) else []
        crossing_rows[name] = hit_class
        filippov[name] = filippov_report(name, GUARD_PARTNER[name])

    primary["_guard_summary"] = guard_summary
    primary["_period"] = period
    returns = return_diagnostic(pack["points"], primary["_fail"], period)
    drop = transient_drop(primary["_beta"], pack["points"], pack["f_hat"], pack["shape"], period)
    hybrid = hybrid_switches(orbit, period, eq)
    for row in hybrid:
        if row.get("hit_state") is not None:
            klass = classify_points(
                np.asarray(row["hit_state"])[None, :],
                primary["_comps"],
                pack["points"],
                primary["spacing"],
            )[0]
            row["hit_class"] = klass["class"]
            row["hit_component_kind"] = klass["component_kind"]
            row["hit_dist_nearest"] = klass["dist_nearest"]
            row["hit_dist_second"] = klass["dist_second"]

    # Verdict from the predeclared rule.
    interior_hits = []
    for name, summary in guard_summary.items():
        if summary["counts"]["interior"] > 0:
            interior_hits.append(name)
    crossing_interior = []
    for name, rows in crossing_rows.items():
        for r in rows:
            if r["class"] == "interior":
                crossing_interior.append({"guard": name, **r})
    geometric_orbit_hits = {name: int(len(crossings[name])) for name in crossings}
    verdict = "counterexample" if (interior_hits or crossing_interior or any(geometric_orbit_hits[n] > 0 and n != "placeholder" for n in geometric_orbit_hits)) else "alignment"
    # Geometric membership in the orbit is already a counterexample to "switches only at interfaces",
    # but only when the hit is not the focus. The focus lies on imm and ang by construction.
    # Require a hit whose distance to the focus exceeds one spacing.
    geometric_off_focus = []
    for name, pts in crossings.items():
        for p in pts:
            if float(np.linalg.norm(p - eq)) > primary["spacing"]:
                geometric_off_focus.append({"guard": name, "T": float(p[0]), "E": float(p[1]), "dist_focus": float(np.linalg.norm(p - eq))})
    if geometric_off_focus or interior_hits or crossing_interior:
        verdict = "counterexample"
    else:
        verdict = "alignment"

    draw(primary, orbit, eq, crossings, hybrid, drop)

    def public_record(rec: dict) -> dict:
        return {
            "n": rec["n"],
            "n_points": rec["n_points"],
            "spacing": rec["spacing"],
            "shape": rec["shape"],
            "support_radius": rec["support_radius"],
            "matrix_symmetry_max": rec["matrix_symmetry_max"],
            "condition_number": rec["condition_number"],
            "min_eigenvalue": rec["min_eigenvalue"],
            "gamma": rec["gamma"],
            "fail_count": rec["fail_count"],
            "fail_fraction": rec["fail_fraction"],
            "orbit_cover_within_2_spacing": rec["orbit_cover_within_2_spacing"],
            "node_residual_inf": rec["node_residual_inf"],
            "flatness": rec["flatness"],
            "components": rec["components"],
            "selected": bool(rec is chosen[rec["n"]]),
        }

    grids_public = {}
    for n, records in by_grid.items():
        grids_public[str(n)] = [public_record(r) for r in records]

    hybrid_public = []
    for row in hybrid:
        hybrid_public.append({k: v for k, v in row.items() if not k.startswith("_")})

    focus_comp = [c for c in primary["components"] if c["kind"] == "focus"]
    orbit_comp = [c for c in primary["components"] if c["kind"] == "orbit"]

    cross_grid = {}
    for n, rec in chosen.items():
        per_guard = {}
        for name, pts in crossings.items():
            rows = classify_points(pts, rec["_comps"], rec["_pack"]["points"], rec["spacing"]) if len(pts) else []
            per_guard[name] = rows
        cross_grid[str(n)] = per_guard

    results = {
        "run_label": RUN_LABEL,
        "deterministic": True,
        "constants": {"a": A_PAR, "b": B_PAR, "c": C_PAR, "d": D_PAR},
        "closed_form": {"T_star": T_STAR, "E_star": E_STAR, "y_det": Y_DET},
        "window": {"T": [T_LO, T_HI], "E": [E_LO, E_HI]},
        "collocation": {
            "kernel": "Wendland psi_4,2",
            "delta2": DELTA2,
            "support_multiplier": SUPPORT_MULT,
            "stencil_each_way": STENCIL_M,
            "stencil_step_in_spacings": STENCIL_STEP,
            "gamma_menu": list(GAMMA_MENU),
            "selection_rule": "Among cuts covering at least 95 percent of the orbit within two spacings, keep the smallest peak-to-peak V. If none qualify, keep the smallest peak-to-peak anyway.",
            "classification": {
                "interior": "nearest defect component within 1 spacing and every other component beyond 2 spacings",
                "interface": "two components of different kind both within 1.5 spacings",
                "transient": "every failing node beyond 2 spacings",
                "margin": "otherwise",
            },
        },
        "equilibria": equilibria,
        "orbit": {"period": period, **orbit_meta},
        "attractiveness": attract,
        "pause_fields_in_window": pause_eq,
        "grids": grids_public,
        "primary_n": N_PRIMARY,
        "primary_gamma": primary["gamma"],
        "binary_iterations_not_reported_function": history,
        "guards": {
            "cyc": {"observable": "T - b", "threshold": B_PAR, "partner": "Q", "summary": guard_summary["cyc"], "orbit_crossings": crossing_rows["cyc"], "filippov": filippov["cyc"]},
            "imm": {"observable": "E - E*", "threshold": E_STAR, "partner": "I", "summary": guard_summary["imm"], "orbit_crossings": crossing_rows["imm"], "filippov": filippov["imm"]},
            "ang": {"observable": "T - T*", "threshold": T_STAR, "partner": "A", "summary": guard_summary["ang"], "orbit_crossings": crossing_rows["ang"], "filippov": filippov["ang"]},
        },
        "geometric_orbit_hits_off_focus": geometric_off_focus,
        "crossing_classes_by_grid": cross_grid,
        "return_diagnostic": returns,
        "transient": {
            "start": drop["start"],
            "V_start": drop["V_start"],
            "V_after_one_period": drop["V_after_one_period"],
            "drop_one_period": drop["drop_one_period"],
            "V_end": drop["V_end"],
            "drop_total": drop["drop_total"],
        },
        "hybrid_switches": hybrid_public,
        "component_counts_primary": {
            "n_connected": len(primary["components"]),
            "n_focus": len(focus_comp),
            "n_orbit": len(orbit_comp),
            "focus_sizes": [c["size"] for c in focus_comp],
            "orbit_sizes": [c["size"] for c in orbit_comp],
        },
        "verdict": verdict,
        "honesty": {
            "collocation_defect_is_not_a_certified_conley_set": True,
            "reason": "The rebuilt orbital derivative along the periodic orbit is not identically non-positive, no epsilon-chain is enumerated, and no interval enclosure is computed.",
            "primary_orbit_fraction_positive_orbital_derivative": primary["flatness"]["fraction_positive"],
        },
    }
    payload = jsonable(results)
    (ROOT / "results.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "verdict": verdict,
        "period": period,
        "primary_gamma": primary["gamma"],
        "fail_fraction": primary["fail_fraction"],
        "fraction_positive": primary["flatness"]["fraction_positive"],
        "V_pp": primary["flatness"]["V_peak_to_peak"],
        "components": payload["component_counts_primary"],
        "guards": {k: payload["guards"][k]["summary"] for k in payload["guards"]},
        "crossings": {k: payload["guards"][k]["orbit_crossings"] for k in payload["guards"]},
        "geometric_off_focus": geometric_off_focus,
        "hybrid": hybrid_public,
        "orbit_box": orbit_meta,
        "eq": equilibria["coexistence"],
        "transient_drop": payload["transient"]["drop_total"],
        "return": returns,
        "filippov_opposite": {k: filippov[k]["fraction_opposite"] for k in filippov},
        "pause": pause_eq,
        "attract": attract,
        "binary": history,
        "coarse": public_record(chosen[N_COARSE])["fail_fraction"],
        "fine": public_record(chosen[N_FINE])["fail_fraction"],
        "cover": {
            "coarse": public_record(chosen[N_COARSE])["orbit_cover_within_2_spacing"],
            "primary": primary["orbit_cover_within_2_spacing"],
            "fine": public_record(chosen[N_FINE])["orbit_cover_within_2_spacing"],
        },
    }, indent=2, default=str))


if __name__ == "__main__":
    main()
