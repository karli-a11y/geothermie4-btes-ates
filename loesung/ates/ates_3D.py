"""ATES 3D MIT regionaler Grundwasserstroemung, i = 0,015.

Was von der Speicherwirkung uebrig bleibt, wenn die natuerliche Stroemung
mitrechnet. Bei kf = 6e-4 m/s und i = 0,015 wandert die Waermefahne rund
490 m im Jahr - was im Sommer eingespeichert wird, steht im Winter einen
halben Kilometer stromab. Deshalb genuegen 1 bis 2 Betriebsjahre: es gibt
nichts, was sich ueber Jahre aufladen koennte.

    python ates_3D.py             # 2 Betriebsjahre, rund 36 min
    python ates_3D.py --years 1   # halbe Rechenzeit
    python ates_3D.py --no-run    # nur Netz und .prj, nicht rechnen

Zum Bearbeiten gibt es genau einen Block: das Dict FALL weiter unten.

Der Lauf legt ergebnisse_3d/ an: ates3d_gw.pvd und die VTU je Zeitschritt fuer
ParaView, csv/ mit Zeitreihe, Monatsbilanz, Jahreskennzahlen und Pruefblatt,
figures/ mit den Abbildungen; 7_draufsicht.png und 8_laengsschnitt.png zeigen
die Drift.

In ParaView die .pvd oeffnen und die Farbskala fuer T fest auf 283,15-333,15 K
stellen, sonst skaliert jeder Zeitschritt neu.
"""
from __future__ import annotations

from pathlib import Path

FALL = {

    "gw_gradient":        0.015,
    "gw_richtung_grad":     0.0,

    "monatsleistung_W": [
        -595_850.0,
        -523_490.0,
        -336_450.0,
          +4_975.0,
        +162_190.0,
        +343_595.0,
        +346_995.0,
        +282_210.0,
         +24_135.0,
        -281_920.0,
        -370_970.0,
        -447_005.0,
    ],

    "T_injektion_C":   60.0,
    "T_aquifer_C":     10.0,
    "betriebsjahre":      2,

    "aquifer": {
        "maechtigkeit_m":                 38.0,
        "kf_m_s":                       6.0e-4,
        "porositaet":                   0.1191,
        "dichte_korn_kg_m3":            2760.0,
        "waermekapazitaet_korn_J_kgK":   793.0,
        "waermeleitfaehigkeit_korn_W_mK": 2.28,
    },

    "deckgestein": {
        "maechtigkeit_m":                60.0,
        "permeabilitaet_m2":          2.1e-16,
        "porositaet":                    0.05,
        "dichte_korn_kg_m3":           2700.0,
        "waermekapazitaet_korn_J_kgK":  900.0,
        "waermeleitfaehigkeit_korn_W_mK": 2.0,
    },

    "filter_kantenlaenge_m": 1.0,

}

#!/usr/bin/env python3
"""ates_report.py — automatischer Prüf- und Auswertebericht nach jedem ATES-Lauf.

Warum es das gibt
-----------------
Ein ATES-Lauf kann numerisch tadellos durchlaufen und trotzdem physikalischen
Unsinn liefern. Zwei Fehler sind dabei besonders heimtückisch, weil man sie den
Ergebnissen nicht ansieht, wenn man die Vergleichszahl nicht kennt:

  A) Der Quellterm der HT-DRUCKgleichung ist eine VOLUMENbilanz [1/s], nicht
     massenbasiert [kg/(m³·s)]. Fehlt die Division durch ρ_f, injiziert das
     Modell den 1000-fachen Volumenstrom.
     Fingerabdruck: v_Darcy einige 100 m/d statt ~1 m/d, Brunnendruck mehrere
     bar statt einiger kPa, fast der ganze Aquifer über 50 °C.

  B) Die Dirichlet-Temperatur am Brunnen darf NUR während der Beladung aktiv
     sein (`DirichletWithinTimeInterval`). Läuft sie auch beim Fördern, ist die
     Entnahmetemperatur vorgegeben statt berechnet und die Energiebilanz kann
     nicht schließen.
     Fingerabdruck: T am Brunnen exakt konstant über Wochen, T_max über T_inj.

Dieses Modul rechnet nach jedem Lauf die Kennzahlen aus, prüft sie gegen
plausible Bänder und schreibt Abbildungen, die die Argumentationskette zeigen
statt nur Kurven. Es funktioniert für das 2D-Radialmodell und die 3D-Modelle
gleichermaßen.

Aufruf
------
Automatisch aus dem Modellskript::

    import ates_report
    ates_report.auto_report(CONFIG, out_dir, curves=curves)

Von Hand für einen schon gerechneten Lauf::

    python ates_report.py            # nimmt CONFIG aus dem Skript im Ordner
    ATES_REPORT=0 python ates_3d.py  # Bericht abschalten

Das Modul importiert KEIN Modellskript (sonst Zirkelimport) — die CONFIG kommt
als Argument. Es fängt jede Ausnahme selbst ab: kein Lauf darf am Bericht
scheitern.
"""

import csv
import math
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

DAY = 86400.0
YEAR = 365.25 * DAY
MONTH = YEAR / 12.0
G = 9.81
MONTHS = ("Jan", "Feb", "Mär", "Apr", "Mai", "Jun",
          "Jul", "Aug", "Sep", "Okt", "Nov", "Dez")

C_BLUE, C_ORANGE, C_AQUA, C_YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
C_GOOD, C_WARN, C_CRIT = "#0d6b4f", "#a1620b", "#b3261e"
C_INK, C_INK2, C_INK3, C_RULE = "#14171a", "#4d565e", "#808b95", "#d8dde2"

@dataclass
class RunInfo:
    out_dir: Path
    prefix: str
    axisym: bool
    snapshots: list = field(default_factory=list)
    complete: bool = True

    @property
    def t_last(self) -> float:
        return self.snapshots[-1][0] if self.snapshots else 0.0

def detect_run(cfg: dict, out_dir=None) -> RunInfo:
    prefix = cfg["output"]["prefix"]
    out = Path(out_dir if out_dir is not None else cfg["output"]["out_dir"])
    snaps = []
    for p in out.glob(f"{prefix}_ts_*.vtu"):
        m = re.search(r"_t_([0-9.]+)\.vtu$", p.name)
        if m:
            snaps.append((float(m.group(1)), p))
    snaps.sort(key=lambda x: x[0])

    axisym = False
    prj = out / f"{prefix}.prj"
    if prj.exists():
        try:
            axisym = 'axially_symmetric="true"' in prj.read_text(
                encoding="ISO-8859-1", errors="replace")
        except OSError:
            pass
    return RunInfo(out_dir=out, prefix=prefix, axisym=axisym, snapshots=snaps)

def well_config(cfg: dict) -> dict:
    w = cfg.get("well") or cfg.get("wells") or {}
    lay = cfg["layers"]
    t_aq = lay["aquifer_thickness_m"]
    off = w.get("screen_top_offset_m", 0.0) + w.get("screen_bottom_offset_m", 0.0)
    h_screen = max(1e-9, t_aq - off)

    if "r_well_m" in w:
        r_eq = w["r_well_m"]
        V = math.pi * r_eq ** 2 * h_screen
    else:
        dx, dy = w.get("screen_dx_m", 1.0), w.get("screen_dy_m", 1.0)
        V = dx * dy * h_screen
        r_eq = math.sqrt(dx * dy / math.pi)

    dom = cfg["domain"]
    R = dom.get("r_max_m") or 0.5 * min(dom.get("size_x_m", 300.0),
                                        dom.get("size_y_m", 300.0))
    return {
        "r_eq_m": r_eq, "h_screen_m": h_screen, "V_well_m3": V,
        "A_screen_m2": 2.0 * math.pi * r_eq * h_screen,
        "R_influence_m": R,
        "production_control": w.get("production_control", "fixed"),
        "max_rate_factor": w.get("max_rate_factor", 1.0),
    }

def cycle_info(cfg: dict, curves: dict | None = None, rate_mult=None) -> dict:
    cyc = cfg["cycles"]
    op = cfg["operation"]
    T_amb = cfg["initial"]["T_K"]
    dT_ref = op["T_hot_K"] - T_amb
    monthly = cyc.get("monthly_power_W")
    cp_f = cfg["fluid"]["cp_J_kgK"]

    if curves:
        mdot_nom = curves.get("mdot_nom_kg_s") or op["mass_flow_rate_kg_s"]
        charge = list(curves.get("charge_intervals") or [])
        power = curves.get("cycle_power")
        t_total = curves.get("t_total")
    else:
        mdot_nom, charge, power, t_total = op["mass_flow_rate_kg_s"], [], None, None
    if monthly and curves is None:
        mdot_nom = max(abs(p) for p in monthly) / (cp_f * dT_ref)

    return {
        "mode": "monthly" if monthly else "phases",
        "monthly_power_W": list(monthly) if monthly else None,
        "mdot_nom_kg_s": mdot_nom, "charge_intervals": charge,
        "cycle_power": power, "t_total": t_total,
        "T_amb_K": T_amb, "dT_ref_K": dT_ref, "cp_f": cp_f,
        "n_cycles": cyc["n_cycles"], "rate_mult": rate_mult,
    }

@dataclass
class Geometry:
    vol: np.ndarray
    rc: np.ndarray
    conn: np.ndarray | None
    mask_aq: np.ndarray
    mask_well: np.ndarray
    r_cell: np.ndarray
    z_cell: np.ndarray
    w_well: np.ndarray
    x_rel: np.ndarray
    z_aq: tuple
    n_points: int
    n_cells: int
    h_typ_well_m: float

def load_geometry(cfg: dict, run: RunInfo):
    import pyvista as pv
    dom = pv.read(run.out_dir / f"{run.prefix}_domain.vtu")
    mid = np.asarray(dom.cell_data["MaterialIDs"])

    mats = {0: cfg["materials"]["aquifer"],
            1: cfg["materials"]["caprock_top"],
            2: cfg["materials"]["caprock_bottom"],
            3: cfg["materials"]["aquifer"],
            4: cfg["materials"]["aquifer"]}
    rho_f, cp_f = cfg["fluid"]["rho_ref_kg_m3"], cfg["fluid"]["cp_J_kgK"]

    def rc_of(i):
        m = mats.get(int(i), cfg["materials"]["aquifer"])
        n = m["porosity"]
        return n * rho_f * cp_f + (1 - n) * m["rho_s_kg_m3"] * m["cp_s_J_kgK"]

    rc = np.array([rc_of(i) for i in mid])

    sizes = dom.compute_cell_sizes(length=False, area=True, volume=True)
    cc = dom.cell_centers().points
    if run.axisym:
        vol = 2.0 * math.pi * cc[:, 0] * np.abs(sizes.cell_data["Area"])
        r_cell, z_cell = cc[:, 0], cc[:, 1]
        x_rel = cc[:, 0]
    else:
        vol = np.abs(sizes.cell_data["Volume"])
        wsel = (mid == 3)
        if wsel.any():
            xw = float(np.average(cc[wsel, 0], weights=np.abs(sizes.cell_data["Volume"])[wsel]))
            yw = float(np.average(cc[wsel, 1], weights=np.abs(sizes.cell_data["Volume"])[wsel]))
        else:
            xw = yw = 0.0
        r_cell = np.hypot(cc[:, 0] - xw, cc[:, 1] - yw)
        x_rel = cc[:, 0] - xw
        z_cell = cc[:, 2]

    keys = list(dom.cells_dict.keys())
    conn = dom.cells_dict[keys[0]] if len(keys) == 1 else None

    z_base = cfg["domain"]["z_base_m"]
    z0 = z_base + cfg["layers"]["caprock_bottom_thickness_m"]
    z1 = z0 + cfg["layers"]["aquifer_thickness_m"]
    mask_well = np.isin(mid, (3,))
    mask_aq = np.isin(mid, (0, 3, 4)) & (z_cell >= z0 - 1e-6) & (z_cell <= z1 + 1e-6)
    if not mask_well.any():
        mask_well = mask_aq & (r_cell < 2.0)
    w_well = vol[mask_well] / vol[mask_well].sum()
    h_typ = float(np.median(vol[mask_well] ** (1 / 3))) if not run.axisym else \
        float(np.median(np.abs(sizes.cell_data["Area"][mask_well]) ** 0.5))
    return dom, Geometry(vol=vol, rc=rc, conn=conn, mask_aq=mask_aq,
                         mask_well=mask_well, r_cell=r_cell, z_cell=z_cell,
                         w_well=w_well, x_rel=x_rel, z_aq=(z0, z1), n_points=dom.n_points,
                         n_cells=dom.n_cells, h_typ_well_m=h_typ)

def _cell_vals(mesh, name, conn):
    if name not in mesh.point_data:
        return None
    if conn is not None:
        return mesh.point_data[name][conn].mean(axis=1)
    return np.asarray(mesh.point_data_to_cell_data()[name])

def build_timeseries(cfg, run, geo, cyc):
    import pyvista as pv
    rho_f = cfg["fluid"]["rho_ref_kg_m3"]
    T_amb = cyc["T_amb_K"]
    z0, z1 = geo.z_aq
    z_mid = 0.5 * (z0 + z1)
    probe_pt = np.array([[0.8, z_mid, 0.0]]) if run.axisym else np.array([[0.8, 0.0, z_mid]])
    probe = pv.PolyData(probe_pt)

    p_ref = None
    out = {k: [] for k in ("t", "T_well", "T_probe", "T_min", "T_max", "dp_well",
                           "E_aq", "E_cr", "r_front", "v_max", "frac_hot",
                           "x_front_down", "x_front_up")}
    for t_s, f in run.snapshots:
        try:
            m = pv.read(f)
        except Exception:
            run.complete = False
            continue
        T = _cell_vals(m, "T", geo.conn)
        if T is None:
            continue
        dE = geo.rc * geo.vol * (T - T_amb)
        out["t"].append(t_s)
        out["T_well"].append(float(np.dot(T[geo.mask_well], geo.w_well)))
        try:
            out["T_probe"].append(float(probe.sample(m)["T"][0]))
        except Exception:
            out["T_probe"].append(np.nan)
        out["T_min"].append(float(T.min()))
        out["T_max"].append(float(T.max()))
        out["E_aq"].append(float(dE[geo.mask_aq].sum()) / 1e9)
        out["E_cr"].append(float(dE[~geo.mask_aq].sum()) / 1e9)
        hot = geo.mask_aq & (T > T_amb + 1.0)
        out["r_front"].append(float(geo.r_cell[hot].max()) if hot.any() else 0.0)
        out["x_front_down"].append(float(geo.x_rel[hot].max()) if hot.any() else 0.0)
        out["x_front_up"].append(float(geo.x_rel[hot].min()) if hot.any() else 0.0)
        out["frac_hot"].append(float(geo.vol[geo.mask_aq & (T > T_amb + 40.0)].sum()
                                     / geo.vol[geo.mask_aq].sum()) * 100.0)
        p = _cell_vals(m, "p", geo.conn)
        if p is None:
            out["dp_well"].append(np.nan)
        else:
            if p_ref is None:
                p_ref = p.copy()
                _gw = cfg.get("regional_gw", {})
                if _gw.get("enable", False) and not run.axisym:
                    _i = float(_gw.get("gradient_m_per_m", 0.0))
                    _a = math.radians(float(_gw.get("direction_deg", 0.0)))
                    _rg = cfg["fluid"]["rho_ref_kg_m3"] * G * _i
                    try:
                        _xw = float(cfg["wells"]["hot_well_xy"][0])
                    except Exception:
                        _xw = 0.0
                    p_ref = p_ref - _rg * math.cos(_a) * (geo.x_rel + _xw)
            d = p - p_ref
            out["dp_well"].append(float(d[geo.mask_well].mean()))
        v = m.point_data.get("darcy_velocity")
        out["v_max"].append(float(np.linalg.norm(np.asarray(v), axis=1).max())
                            if v is not None else np.nan)

    ts = {k: np.asarray(v, dtype=float) for k, v in out.items()}
    ts["days"] = ts["t"] / DAY

    if cyc["cycle_power"] is not None:
        tg, g = cyc["cycle_power"]
        ts["mdot"] = np.interp(ts["t"], tg, g) * cyc["mdot_nom_kg_s"]
    else:
        ts["mdot"] = np.zeros_like(ts["t"])
    T_inj_t = np.full_like(ts["t"], T_amb)
    for t0, t1, Ti in cyc["charge_intervals"]:
        T_inj_t[(ts["t"] >= t0) & (ts["t"] <= t1)] = Ti
    ts["T_inj"] = T_inj_t
    cp = cyc["cp_f"]
    ts["P"] = np.where(ts["mdot"] > 0,
                       ts["mdot"] * cp * (T_inj_t - T_amb),
                       ts["mdot"] * cp * (ts["T_well"] - T_amb))

    def cumtrap(y):
        if len(ts["t"]) < 2:
            return np.zeros_like(y)
        return np.concatenate([[0.0], np.cumsum(np.diff(ts["t"]) * 0.5 * (y[1:] + y[:-1]))]) / 1e9

    ts["E_in"] = cumtrap(np.clip(ts["P"], 0, None))
    ts["E_out"] = cumtrap(np.clip(-ts["P"], 0, None))

    if len(ts["dp_well"]) and cfg.get("regional_gw", {}).get("enable", False)             and not run.axisym:
        ts["dp_well"][0] = np.nan
    return ts

def energy_metrics(ts, cyc):
    rows = []
    E_tot = ts["E_aq"] + ts["E_cr"]
    monthly = cyc["monthly_power_W"]
    e_dem = (sum(-p for p in monthly if p < 0) * MONTH / 1e9) if monthly else None
    for y in range(int(cyc["n_cycles"])):
        a, b = y * YEAR, (y + 1) * YEAR
        if ts["t"][-1] < b - 0.02 * YEAR:
            break
        e_in = np.interp(b, ts["t"], ts["E_in"]) - np.interp(a, ts["t"], ts["E_in"])
        e_out = np.interp(b, ts["t"], ts["E_out"]) - np.interp(a, ts["t"], ts["E_out"])
        e_field = sum(np.interp(min(t1, ts["t"][-1]), ts["t"], E_tot)
                      - np.interp(t0, ts["t"], E_tot)
                      for t0, t1, _ in cyc["charge_intervals"] if a <= t0 < b)
        sel = (ts["t"] >= a) & (ts["t"] < b)
        prod = sel & (ts["mdot"] < 0)
        w = -ts["mdot"][prod]
        rows.append({
            "jahr": y + 1,
            "E_ein_GJ": e_in, "E_ein_feld_GJ": e_field, "E_aus_GJ": e_out,
            "E_bedarf_GJ": e_dem,
            "eta_pct": 100 * e_out / e_in if e_in > 0 else np.nan,
            "eta_feld_pct": 100 * e_out / e_field if e_field > 0 else np.nan,
            "deckung_pct": (100 * e_out / e_dem) if e_dem else np.nan,
            "T_foerder_C": (float(np.average(ts["T_well"][prod], weights=w)) - 273.15
                            if w.sum() > 0 else np.nan),
            "E_deckgestein_GJ": float(ts["E_cr"][sel].max()) if sel.any() else np.nan,
            "r_front_m": float(ts["r_front"][sel].max()) if sel.any() else np.nan,
        })
    return rows, e_dem

def monthly_balance(ts, cyc, year_idx):
    monthly = cyc["monthly_power_W"]
    if not monthly:
        return []
    rows = []
    for mo in range(12):
        a = year_idx * YEAR + mo * MONTH
        b = a + MONTH
        sel = (ts["t"] >= a) & (ts["t"] < b)
        if not sel.any():
            continue
        prod = sel & (ts["mdot"] < 0)
        e_del = float(np.trapezoid(np.where(prod, np.clip(-ts["P"], 0, None), 0.0),
                                   ts["t"])) / 1e9
        e_dem = max(0.0, -monthly[mo]) * MONTH / 1e9
        w = -ts["mdot"][prod]
        rows.append({
            "monat": MONTHS[mo], "P_soll_kW": monthly[mo] / 1e3,
            "T_well_C": (float(np.average(ts["T_well"][prod], weights=w)) - 273.15
                         if w.sum() > 0 else np.nan),
            "E_gefordert_GJ": e_dem, "E_geliefert_GJ": e_del,
            "deckung_pct": (100 * e_del / e_dem) if e_dem > 0 else np.nan,
            "laden": monthly[mo] > 0,
        })
    return rows

def log_stats(run):
    """Liest die OGS-Konsolenausgabe aus dem Ergebnisverzeichnis."""
    out = {"accepted": None, "rejected": None, "aborted": False, "errors": 0}
    for name in ("driver.log", "run.log", "ogs.log"):
        p = run.out_dir / name
        if not p.exists():
            continue
        try:
            txt = p.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        m = re.findall(r"accepted steps are (\d+), and the rejected steps are (\d+)", txt)
        if m:
            out["accepted"], out["rejected"] = int(m[-1][0]), int(m[-1][1])
        out["errors"] = txt.count("error: ")
        if "terminated with error" in txt or "cannot reduce the time step" in txt:
            out["aborted"] = True
        break
    return out

def checks(cfg, run, geo, ts, cyc, rows):
    """-> Liste von (Name, Wert, Einheit, Status, Diagnose, ok_band, fehler_marke)"""
    T_amb = cyc["T_amb_K"]
    T_inj = cfg["operation"]["T_hot_K"]
    out = []

    def add(name, val, unit, band, status, diag, marke=None):
        out.append(dict(name=name, val=val, unit=unit, band=band,
                        status=status, diag=diag, marke=marke))

    _wc = well_config(cfg)
    v_soll = cyc["mdot_nom_kg_s"] / (cfg["fluid"]["rho_ref_kg_m3"] * _wc["A_screen_m2"])
    v = np.nanmax(ts["v_max"]) if len(ts["v_max"]) else np.nan
    ratio = v / max(v_soll, 1e-30)
    st = "OK" if ratio <= 10 else ("WARNUNG" if ratio <= 100 else "FEHLER")
    add("v_Darcy max / v_erwartet", ratio, "-", (0.0, 10.0), st,
        f"Plausibel ({v*DAY:.2f} m/d gegen erwartete {v_soll*DAY:.2f} m/d)."
        if st == "OK" else
        f"{v*DAY:.0f} m/d gegen erwartete {v_soll*DAY:.2f} m/d - Faktor "
        f"{ratio:.0f}. Bei ~1000 ist der Quellterm der Druckgleichung nicht "
        "durch rho_f geteilt (muss 1/s sein, nicht kg/(m³·s)).", 1000.0)

    K_h = (cfg["materials"]["aquifer"]["permeability_m2"]
           * cfg["fluid"]["rho_ref_kg_m3"] * G / cfg["fluid"]["viscosity_Pa_s"])
    b_aq = cfg["layers"]["aquifer_thickness_m"]
    thiem = (math.log(max(_wc["R_influence_m"] / _wc["r_eq_m"], 1.1))
             / (2 * math.pi * K_h * b_aq))
    dp_soll = (cyc["mdot_nom_kg_s"] / cfg["fluid"]["rho_ref_kg_m3"]
               * thiem * cfg["fluid"]["rho_ref_kg_m3"] * G)
    _dpw = ts["dp_well"][1:] if len(ts["dp_well"]) > 1 else ts["dp_well"]
    dp = np.nanmax(np.abs(_dpw)) if len(_dpw) else np.nan
    q = dp / max(dp_soll, 1e-30)
    st = "OK" if 0.2 <= q <= 4 else ("WARNUNG" if q <= 20 else "FEHLER")
    add("Brunnendruck / Thiem-Vorhersage", q, "-", (0.2, 4.0), st,
        f"Feld ({dp/1e3:.2f} kPa) und Handformel ({dp_soll/1e3:.2f} kPa) "
        "passen zusammen." if st == "OK" else
        f"Feld {dp/1e3:.1f} kPa gegen Handformel {dp_soll/1e3:.1f} kPa - Faktor "
        f"{q:.0f}. Bei ~1000 ist der Quellterm nicht durch rho_f geteilt.", 1000.0)

    fh = np.nanmax(ts["frac_hot"]) if len(ts["frac_hot"]) else np.nan
    st = "OK" if fh <= 15 else ("WARNUNG" if fh <= 40 else "FEHLER")
    add("Aquifervolumen ueber 50 GradC", fh, "%", (0.0, 15.0), st,
        "Fahne bleibt lokal." if st == "OK" else
        "Fast das ganze Gebiet ist heiß — es wird viel zu viel Energie "
        "eingetragen (rho_f) oder die Domäne ist zu klein.", 93.0)

    over = (np.nanmax(ts["T_max"]) - T_inj) if len(ts["T_max"]) else np.nan
    _flow = bool(cfg.get("regional_gw", {}).get("enable", False)) and not run.axisym
    _lift = T_inj - T_amb
    _warn = max(1.0, 0.05 * _lift) if _flow else 1.0
    st = "OK" if over <= 0.05 else ("WARNUNG" if over <= _warn else "FEHLER")
    add("T_max - T_inj", over, "K", (-50.0, 0.05), st,
        "Kein Überschwinger." if st == "OK" else
        (f"Überschwinger am geklemmten Filterkörper im Durchstrom "
         f"({100*over/max(_lift, 1e-9):.1f} % des Temperaturhubs). Bekannter "
         f"Galerkin-Effekt, energetisch belanglos. Erst jenseits von "
         f"{_warn:.1f} K deutet das auf eine Dauerklemme hin."
         if _flow and st == "WARNUNG" else
         "T übersteigt die Injektionstemperatur — es gibt keine Wärmequelle "
         "über T_inj. Advektiver Überschwinger (Netz/Dispersivität) oder "
         "Dauerklemme am Brunnen."), 4.85)

    prod = ts["mdot"] < 0
    sd = float(np.std(ts["T_well"][prod])) if prod.sum() > 3 else np.nan
    st = "OK" if (np.isnan(sd) or sd > 0.5) else "FEHLER"
    add("Streuung T_Brunnen (Förderung)", sd, "K", (0.5, 50.0), st,
        "Fördertemperatur wird berechnet." if st == "OK" else
        "T am Brunnen ist beim Fördern praktisch konstant — die Dirichlet-BC "
        "läuft auch in der Förderphase. `DirichletWithinTimeInterval` nur "
        "über die Beladungsintervalle setzen.", 0.0)

    ls = log_stats(run)
    if ls["accepted"] is not None:
        rej = ls["rejected"]
        acc = ls["accepted"]
        st = ("FEHLER" if ls["aborted"] else
              ("OK" if rej == 0 else ("WARNUNG" if rej <= 0.02 * max(acc, 1) else "FEHLER")))
        add("verworfene Zeitschritte", rej, "-", (0.0, 0.0), st,
            f"{acc} Schritte, keiner verworfen." if st == "OK" else
            (f"Der Lauf ist ABGEBROCHEN ({acc} Schritte, {rej} verworfen). Die "
             "Auswertung unten deckt nur den gerechneten Teil ab. Meist reicht "
             "der Picard-Iteration das Budget nicht: `solver.nonlinear_iter` "
             "erhoehen (das Zeitschritt-Verkleinern hilft dabei NICHT)."
             if ls["aborted"] else
             f"{rej} von {acc} Zeitschritten mussten wiederholt werden - der "
             "Lauf ist durch, aber die Nichtlinearitaet ist grenzwertig. "
             "`solver.nonlinear_iter` erhoehen."))

    und = (T_amb - np.nanmin(ts["T_min"])) if len(ts["T_min"]) else np.nan
    st = "OK" if und <= 1.0 else ("WARNUNG" if und <= 8.0 else "FEHLER")
    add("T_min unter T_amb", und, "K", (0.0, 1.0), st,
        "Keine Unterschwinger." if st == "OK" else
        "Unterschwinger an der Wärmefront. Bis ~8 K ist das der bekannte "
        "Effekt der konsistenten Massenmatrix (tritt auf, solange "
        "dt < h²/(6a)) und energetisch belanglos — kein Modellfehler.", 7.5)

    lam = cfg["materials"]["caprock_top"]
    n = lam["porosity"]
    a_cr = ((n * cfg["fluid"]["lambda_W_mK"] + (1 - n) * lam["lambda_s_W_mK"])
            / (n * cfg["fluid"]["rho_ref_kg_m3"] * cfg["fluid"]["cp_J_kgK"]
               + (1 - n) * lam["rho_s_kg_m3"] * lam["cp_s_J_kgK"]))
    t_run = max(ts["t"][-1], 1.0) if len(ts["t"]) else 1.0
    front = 2.0 * math.sqrt(a_cr * t_run)
    d_cr = min(cfg["layers"]["caprock_top_thickness_m"],
               cfg["layers"]["caprock_bottom_thickness_m"])
    st = ("OK" if d_cr >= 1.5 * front
          else ("WARNUNG" if d_cr >= front else "FEHLER"))
    add("Deckgestein / Leitfront 2*sqrt(a*t)", d_cr / max(front, 1e-9), "-",
        (1.5, 10.0), st,
        f"Deckgestein {d_cr:.0f} m > Front {front:.0f} m." if st == "OK" else
        f"Die Wärmeleitfront ({front:.0f} m) erreicht den auf T_amb fixierten "
        f"Rand ({d_cr:.0f} m) — `layers.caprock_*_thickness_m` erhöhen, sonst "
        "werden die Verluste künstlich groß.")

    gw_on = bool(cfg.get("regional_gw", {}).get("enable", False)) and not run.axisym
    if gw_on:
        _dn = float(np.nanmax(ts["x_front_down"])) if len(ts.get("x_front_down", [])) else np.nan
        try:
            _xw = float(cfg["wells"]["hot_well_xy"][0])
        except Exception:
            _xw = 0.0
        _av = 0.5 * cfg["domain"].get("size_x_m", 400.0) - _xw
        st = "OK" if _dn <= 0.85 * _av else "WARNUNG"
        add("Fahne stromab / Platz stromab", _dn / max(_av, 1e-9), "-", (0.0, 0.85), st,
            f"Fahne {_dn:.0f} m stromab bei {_av:.0f} m Platz - die Drift bleibt sichtbar."
            if st == "OK" else
            f"Die Fahne ({_dn:.0f} m) erreicht den Ausstroemrand ({_av:.0f} m). Fuer die "
            "Aussage 'sie laeuft weg' reicht das; wer die Fahnenform bis zum Ende sehen "
            "will, braucht mehr `domain.size_x_m` oder muss den Brunnen weiter stromauf setzen.")
    else:
        rf = np.nanmax(ts["r_front"]) if len(ts["r_front"]) else np.nan
        Rd = well_config(cfg)["R_influence_m"]
        st = "OK" if rf <= 0.6 * Rd else ("WARNUNG" if rf <= 0.9 * Rd else "FEHLER")
        add("Fahnenreichweite / Modellrand", rf / max(Rd, 1e-9), "-", (0.0, 0.6), st,
            f"Fahne {rf:.0f} m, Rand {Rd:.0f} m." if st == "OK" else
            f"Die 1-K-Front ({rf:.0f} m) kommt dem Modellrand ({Rd:.0f} m) zu nahe "
            "— Domäne vergrößern (`domain.r_max_m` bzw. `size_x_m/size_y_m`).")

    if len(rows) >= 4:
        e = [r["eta_pct"] for r in rows]
        d = (e[-1] + e[-2]) / 2 - (e[-3] + e[-4]) / 2
        st = "OK" if abs(d) <= 2.0 else "WARNUNG"
        add("eta-Trend (2-Jahres-Mittel)", d, "pp", (-2.0, 2.0), st,
            "Eingeschwungen." if st == "OK" else
            "Noch nicht eingeschwungen — `cycles.n_cycles` erhöhen, sonst ist "
            "der Rückgewinnungsgrad nicht aussagekräftig.")

    wc = well_config(cfg)
    mdot_max = cyc["mdot_nom_kg_s"] * max(1.0, wc["max_rate_factor"])
    v_entry = mdot_max / (cfg["fluid"]["rho_ref_kg_m3"] * wc["A_screen_m2"])
    st = "OK" if v_entry <= 0.01 else ("WARNUNG" if v_entry <= 0.03 else "FEHLER")
    add("Filtereintrittsgeschw.", v_entry, "m/s", (0.0, 0.03), st,
        "Weit unter dem Richtwert 0,03 m/s." if st == "OK" else
        "Über dem Richtwert 0,03 m/s (Sandeintrag, Verkrustung) — Filter "
        "verlängern/aufweiten oder Rate senken.")

    if rows:
        r = rows[-1]
        rest = r["E_ein_GJ"] - r["E_aus_GJ"] - (r["E_ein_feld_GJ"] - r["E_aus_GJ"])
        rel = 100 * (r["E_ein_feld_GJ"] - r["E_ein_GJ"]) / max(r["E_ein_GJ"], 1e-9)
        if gw_on:
            st, lim = "OK", (0.0, 400.0)
            txt = ("Waermetauschereffekt am Grundwasserstrom - so erwartet. "
                   "Deshalb ist eta im Stroemungsfall wenig aussagekraeftig; "
                   "belastbar ist dort die Foerdertemperatur.")
        else:
            st = "OK" if rel <= 30 else ("WARNUNG" if rel <= 60 else "FEHLER")
            lim = (0.0, 30.0)
            txt = ("Der Brunnen traegt etwas mehr ein als mdot*c_p*dT - "
                   "normal, skaliert mit dem geklemmten Filtervolumen."
                   if st == "OK" else
                   "Der Dirichlet-Block traegt sehr viel mehr ein als "
                   "mdot*c_p*dT - Filtervolumen verkleinern "
                   "(`screen_dx_m`/`screen_dy_m`).")
        add("Dirichlet-Mehreintrag", rel, "%", lim, st, txt)
    return out

def _style(ax, title=None, sub=None, xl=None, yl=None):
    ax.grid(alpha=.25, lw=.7, color=C_RULE)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(C_RULE)
    ax.tick_params(colors=C_INK2, labelsize=9)
    if title:
        ax.set_title(title, color=C_INK, fontsize=12, fontweight="bold",
                     loc="left", pad=26)
    if sub:
        ax.text(0, 1.02, sub, transform=ax.transAxes, color=C_INK3,
                fontsize=9.5, va="bottom")
    if xl:
        ax.set_xlabel(xl, color=C_INK2, fontsize=10)
    if yl:
        ax.set_ylabel(yl, color=C_INK2, fontsize=10)

def fig_pruefblatt(chks, fig_dir, plt):
    """Ampel-Panel: liegt jede Kennzahl im plausiblen Band?"""
    ch = [c for c in chks if np.isfinite(c["val"])]
    if not ch:
        return None
    n = len(ch)
    fig, ax = plt.subplots(figsize=(12.6, 0.62 * n + 2.4))
    col = {"OK": C_GOOD, "WARNUNG": C_WARN, "FEHLER": C_CRIT}
    for i, c in enumerate(reversed(ch)):
        y = i
        lo, hi = c["band"]
        v = c["val"]
        span = max(hi - lo, 1e-12)
        xn = (v - lo) / span
        ax.axhspan(y - .34, y + .34, xmin=0, xmax=1, color="#f4f6f7", lw=0)
        ax.barh(y, 1.0, left=0, height=.34, color=C_GOOD, alpha=.16, lw=0)
        clipped = not (-0.28 <= xn <= 1.28)
        xp = min(max(xn, -0.26), 1.26)
        ax.plot([xp], [y], ">" if xn > 1.28 else ("<" if xn < -0.28 else "o"),
                ms=10 if clipped else 9, color=col[c["status"]], zorder=5)
        ax.text(-0.32, y, c["name"], ha="right", va="center",
                fontsize=9.5, color=C_INK)
        ax.text(1.36, y, f"{v:,.3g} {c['unit']}".replace(",", " "),
                ha="left", va="center", fontsize=9.5, color=C_INK,
                fontfamily="monospace")
        if c["marke"] is not None:
            ax.text(1.95, y, f"Fehlerbild: {c['marke']:g}", ha="left",
                    va="center", fontsize=8.5, color=C_CRIT, alpha=.75)
        ax.text(2.52, y, c["status"], ha="left", va="center", fontsize=9,
                color=col[c["status"]], fontweight="bold")
    ax.set_xlim(-0.35, 3.05)
    ax.set_ylim(-.8, n - .2)
    ax.set_yticks([])
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["untere\nOK-Grenze", "obere\nOK-Grenze"], fontsize=8.5,
                       color=C_INK3)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(C_RULE)
    ax.axvline(0, color=C_GOOD, lw=1, alpha=.5)
    ax.axvline(1, color=C_GOOD, lw=1, alpha=.5)
    fig.suptitle("Prüfblatt: ist dieser Lauf gesund?", x=.012, ha="left",
                 fontsize=13, fontweight="bold", color=C_INK)
    fig.text(.012, .945, "grünes Band = plausibel · Punkt = dieser Lauf · "
             "Pfeil = liegt außerhalb der Skala · rechts das Fehlerbild zum Vergleich",
             fontsize=9, color=C_INK3, va="top")
    bad = [c for c in ch if c["status"] != "OK"]
    txt = ("Alle Kennzahlen im plausiblen Band."
           if not bad else "\n".join(f"• {c['name']}: {c['diag']}" for c in bad[:4]))
    fig.text(.012, .012, txt, fontsize=8.6, color=(C_INK if not bad else C_CRIT),
             va="bottom", wrap=True)
    fig.tight_layout(rect=(0, .11 if bad else .04, 1, .93))
    p = fig_dir / "0_pruefblatt.png"
    fig.savefig(p, dpi=140); plt.close(fig)
    return p

def fig_brunnentemperatur(cfg, ts, cyc, fig_dir, plt):
    T_amb, T_inj = cyc["T_amb_K"] - 273.15, cfg["operation"]["T_hot_K"] - 273.15
    yr = ts["days"] / 365.25
    fig, ax = plt.subplots(figsize=(11.5, 4.4))
    ax.axhline(T_inj, color=C_INK3, lw=1.2, ls="--")
    ax.axhline(T_amb, color=C_INK3, lw=1.2, ls="--")
    ax.text(yr[-1], T_inj, " T_inj", va="center", fontsize=9, color=C_INK3)
    ax.text(yr[-1], T_amb, " T_amb", va="center", fontsize=9, color=C_INK3)
    for t0, t1, _ in cyc["charge_intervals"]:
        ax.axvspan(t0 / YEAR, t1 / YEAR, color=C_ORANGE, alpha=.07, lw=0)
    if np.isfinite(ts["T_probe"]).any():
        ax.plot(yr, ts["T_probe"] - 273.15, lw=1.0, color=C_INK3, ls=":",
                label="Punktsonde halbe Aquiferhöhe (liest zu kalt)")
    ax.plot(yr, ts["T_well"] - 273.15, lw=1.5, color=C_ORANGE,
            label="T am Brunnen (Mittel über die Filterstrecke)")
    ax.plot(yr, ts["T_min"] - 273.15, lw=.8, color=C_BLUE, alpha=.7,
            label="T_min im Gebiet")
    prod = ts["mdot"] < 0
    if prod.sum() > 3:
        sd = float(np.std(ts["T_well"][prod]))
        d = float(np.nanmean(ts["T_well"] - ts["T_probe"]))
        note = f"Streuung von T beim Fördern: {sd:.2f} K"
        if sd < 0.05:
            note += "  ← praktisch konstant: Dirichlet läuft auch beim Fördern!"
        if np.isfinite(d):
            note += f"\nFiltermittel − Punktsonde: {d:+.2f} K (Auftrieb)"
        ax.text(.012, .04, note, transform=ax.transAxes, fontsize=9,
                color=(C_CRIT if sd < 0.05 else C_INK2), va="bottom",
                bbox=dict(fc="white", ec=C_RULE, lw=.8, pad=4))
    _style(ax, "Brunnentemperatur: Vorgabe oder Ergebnis?",
           "Beim Beladen ist T vorgegeben (rot hinterlegt). Beim Fördern muss sie "
           "sich einstellen — eine flache Linie wäre ein Modellfehler.",
           "Betriebsjahr", "T [°C]")
    ax.legend(fontsize=8.5, frameon=False, loc="upper right")
    fig.tight_layout()
    p = fig_dir / "1_brunnentemperatur.png"
    fig.savefig(p, dpi=140); plt.close(fig)
    return p

def fig_deckungsgrad(rows, cyc, fig_dir, plt):
    monthly = cyc["monthly_power_W"]
    if not rows:
        return None
    j = [r["jahr"] for r in rows]
    fig, ax = plt.subplots(figsize=(10.5, 4.6))
    cap = None
    if monthly:
        ein = sum(p for p in monthly if p > 0)
        aus = sum(-p for p in monthly if p < 0)
        cap = 100 * ein / aus if aus > 0 else None
    if cap:
        ax.axhline(cap, color=C_CRIT, lw=1.6)
        ax.annotate(f"Obergrenze des Deckungsgrads: {cap:.1f} %\n"
                    "(mehr Wärme als eingespeichert kann nicht heraus)",
                    xy=(j[len(j)//3], cap), xytext=(0, -34),
                    textcoords="offset points", fontsize=9, color=C_CRIT,
                    arrowprops=dict(arrowstyle="->", color=C_CRIT, lw=1.1))
    ax.axhspan(50, 80, facecolor=C_AQUA, alpha=.05, lw=0)
    ax.text(j[-1], 78, " η realer\n ATES", fontsize=8.5, color=C_INK3, va="top")
    ax.plot(j, [r["eta_pct"] for r in rows], "s-", ms=4, lw=1.4, color=C_ORANGE,
            label="η = E_aus/E_ein")
    ax.plot(j, [r["eta_feld_pct"] for r in rows], "^-", ms=4, lw=1.2,
            color=C_YELLOW, label="η feldbasiert")
    if monthly:
        d = [r["deckung_pct"] for r in rows]
        ax.plot(j, d, "o-", ms=5, lw=2.0, color=C_BLUE, label="Deckungsgrad")
        ax.fill_between(j, 0, d, color=C_BLUE, alpha=.10)
        ax.annotate(f"{d[-1]:.1f} %", (j[-1], d[-1]), textcoords="offset points",
                    xytext=(6, -2), fontsize=10, color=C_BLUE, fontweight="bold")
        if cap:
            ax.text(.012, .96, f"Deckungsgrad = η · {cap:.1f} %",
                    transform=ax.transAxes, fontsize=10, color=C_INK, va="top",
                    bbox=dict(fc="white", ec=C_RULE, lw=.8, pad=4))
    ax.set_ylim(0, 100)
    _style(ax, "Trägt der Speicher die Last?",
           "η sagt, wie gut der Speicher arbeitet. Der Deckungsgrad sagt, ob es "
           "reicht — und der ist durch das Lastprofil nach oben begrenzt."
           if monthly else
           "4-Phasen-Modus: ohne Lastprofil ist kein Deckungsgrad definiert.",
           "Betriebsjahr", "[%]")
    ax.legend(fontsize=9, frameon=False, loc="lower right")
    fig.tight_layout()
    p = fig_dir / "2_deckungsgrad.png"
    fig.savefig(p, dpi=140); plt.close(fig)
    return p

def fig_monatsbilanz(mon, cyc, fig_dir, plt):
    if not mon:
        return None
    x = np.arange(len(mon))
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(10.5, 5.6), sharex=True,
                                  gridspec_kw=dict(height_ratios=[2, 1]))
    dem = [m["E_gefordert_GJ"] for m in mon]
    dele = [m["E_geliefert_GJ"] for m in mon]
    ax.bar(x, dem, .62, facecolor="none", edgecolor=C_INK3, lw=1.1,
           label="gefordert (Lastprofil)")
    ax.bar(x, dele, .62, color=C_BLUE, label="geliefert (Simulation)")
    for i, m in enumerate(mon):
        if np.isfinite(m["deckung_pct"]):
            ax.text(i, max(dem[i], dele[i]) * 1.02, f"{m['deckung_pct']:.0f} %",
                    ha="center", fontsize=8.5, color=C_INK2)
        elif m["laden"]:
            ax.text(i, 2, "laden", ha="center", fontsize=8, color=C_ORANGE,
                    rotation=90, va="bottom")
    _style(ax, "Monatsbilanz im letzten Betriebsjahr",
           "Wo fällt der Speicher aus — und warum?", None, "Energie [GJ]")
    ax.legend(fontsize=9, frameon=False)
    Tw = [m["T_well_C"] for m in mon]
    ax2.plot(x, Tw, "o-", ms=5, lw=1.3, color=C_ORANGE)
    ax2.axhline(cyc["T_amb_K"] - 273.15, color=C_INK3, ls="--", lw=1.1)
    ax2.text(len(mon) - .4, cyc["T_amb_K"] - 273.15, " T_amb", fontsize=8.5,
             color=C_INK3, va="center")
    worst = min((m for m in mon if np.isfinite(m["T_well_C"])),
                key=lambda m: m["T_well_C"], default=None)
    if worst:
        i = [m["monat"] for m in mon].index(worst["monat"])
        ax2.annotate(f"{worst['T_well_C']:.1f} °C → nur "
                     f"{worst['T_well_C'] - (cyc['T_amb_K']-273.15):.1f} K nutzbar",
                     (i, worst["T_well_C"]), textcoords="offset points",
                     xytext=(8, 6), fontsize=9, color=C_CRIT)
    _style(ax2, None, None, None, "T am Brunnen [°C]")
    ax2.set_xticks(x)
    ax2.set_xticklabels([m["monat"] for m in mon])
    fig.tight_layout()
    p = fig_dir / "3_monatsbilanz.png"
    fig.savefig(p, dpi=140); plt.close(fig)
    return p

def fig_machbarkeit(cfg, ts, cyc, mon, fig_dir, plt):
    """Machbarkeitskette: T_inj -> mdot -> Eintrittsgeschw. -> Absenkung."""
    wc = well_config(cfg)
    monthly = cyc["monthly_power_W"]
    if not monthly:
        return None
    rho_f = cfg["fluid"]["rho_ref_kg_m3"]
    mu = cfg["fluid"]["viscosity_Pa_s"]
    k = cfg["materials"]["aquifer"]["permeability_m2"]
    K = k * rho_f * G / mu
    b_aq = cfg["layers"]["aquifer_thickness_m"]
    thiem = math.log(max(wc["R_influence_m"] / wc["r_eq_m"], 1.1)) / (2 * math.pi * K * b_aq)
    cp, dT = cyc["cp_f"], cyc["dT_ref_K"]
    V_LIM = 0.03

    inj = [(MONTHS[i], p / (cp * dT)) for i, p in enumerate(monthly) if p > 0]
    need = []
    for m in mon:
        if not np.isfinite(m["T_well_C"]):
            continue
        d = m["T_well_C"] - (cyc["T_amb_K"] - 273.15)
        if d > .05:
            need.append((m["monat"], abs(m["P_soll_kW"]) * 1e3 / (cp * d)))
    if not (inj or need):
        return None
    lbl = [a for a, _ in inj] + [a for a, _ in need]
    val = np.array([v for _, v in inj] + [v for _, v in need])
    colr = [C_ORANGE] * len(inj) + [C_BLUE] * len(need)
    x = np.arange(len(val))

    def box(ax, lines, color, right=False):
        ax.text(.97 if right else .03, .97, "\n".join(lines),
                transform=ax.transAxes, fontsize=8.8, color=color, va="top",
                ha="right" if right else "left",
                bbox=dict(fc="white", ec=C_RULE, lw=.7, pad=3.5))

    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.9))

    ax = axes[0]
    ax.bar(x, val, .64, color=colr)
    lo, hi = 50 * rho_f / 3600, 250 * rho_f / 3600
    ax.set_ylim(0, val.max() * 1.6)
    for _y, _t in ((lo, "50 m³/h"), (hi, "250 m³/h")):
        if _y < val.max() * 1.5:
            ax.axhline(_y, color=C_AQUA, lw=1.1, ls="--")
            ax.text(len(val) - .5, _y, _t + " ", ha="right", va="bottom",
                    fontsize=8.5, color=C_AQUA)
    box(ax, ["Injektion (orange) folgt aus der ANNAHME T_inj:",
             f"ṁ = P/(c_p·ΔT),  ΔT = {dT:.0f} K",
             "Förderung (blau) = was nötig wäre, um die Last zu decken",
             "gestrichelt: reale ATES-Brunnen 50–250 m³/h"], C_INK2)
    ax.set_title("① Massenstrom", loc="left", fontsize=11.5,
                 fontweight="bold", color=C_INK, pad=10)
    ax.set_ylabel("ṁ [kg/s]", color=C_INK2, fontsize=10)

    ax = axes[1]
    ve = val / (rho_f * wc["A_screen_m2"])
    worst = 100 * ve.max() / V_LIM
    ax.bar(x, 100 * ve / V_LIM, .64, color=colr)
    ax.axhline(100, color=C_CRIT, lw=1.5)
    ax.set_ylim(0, max(115, worst * 1.3))
    ax.text(-.4, 100, "Richtwert 0,03 m/s ", ha="left", fontsize=8.8,
            color=C_CRIT, va="top")
    box(ax, ["v = Q/(2π·r_w·h_Filter)",
             f"Maximum {ve.max():.2e} m/s = {worst:.1f} % des Richtwerts",
             f"Reserve: Faktor {100/max(worst,1e-9):.0f}"],
        C_GOOD if worst < 100 else C_CRIT, right=True)
    ax.set_title("② Beanspruchung des Filters", loc="left", fontsize=11.5,
                 fontweight="bold", color=C_INK, pad=10)
    ax.set_ylabel("in % des Richtwerts", color=C_INK2, fontsize=10)

    ax = axes[2]
    dh = val / rho_f * thiem
    ax.bar(x, dh, .64, color=colr)
    dp_meas = (np.nanmax(np.abs(ts["dp_well"])) / (rho_f * G)
               if len(ts["dp_well"]) else np.nan)
    ax.set_ylim(0, dh.max() * 1.5)
    if np.isfinite(dp_meas):
        ax.axhline(dp_meas, color=C_INK3, ls=":", lw=1.4)
        ax.text(-.4, dp_meas, f"aus dem p-Feld: {dp_meas:.2f} m",
                ha="left", fontsize=8.8, color=C_INK3, va="bottom")
    box(ax, ["Thiem:  dh = Q·ln(R/r_w)/(2π·K·b)",
             f"K = {K:.2e} m/s,  b = {b_aq:.0f} m,  r_w = {wc['r_eq_m']:.2f} m",
             f"Maximum {dh.max():.2f} m",
             "Handformel und p-Feld müssen zusammenpassen —",
             "sonst stimmt der Quellterm nicht."], C_INK2)
    ax.set_title("③ Absenkung / Aufhöhung", loc="left", fontsize=11.5,
                 fontweight="bold", color=C_INK, pad=10)
    ax.set_ylabel("dh [m]", color=C_INK2, fontsize=10)

    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(lbl, fontsize=8.5, rotation=45)
        ax.grid(alpha=.25, lw=.7, color=C_RULE, axis="y")
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(C_RULE)
        ax.tick_params(colors=C_INK2, labelsize=9)

    fig.suptitle("Ist die Anlage baubar? Die Kette T_inj → ṁ → Filter → Druck",
                 x=.008, ha="left", fontsize=13, fontweight="bold", color=C_INK)
    fig.text(.008, .928, "Jede Stufe folgt aus der vorigen. Reißt eine Grenze, ist "
             "die Anlage so nicht baubar. Ein Druck um Faktor 1000 zu hoch heißt: "
             "Quellterm nicht durch ρ_f geteilt.",
             fontsize=9, color=C_INK3, va="top")
    fig.tight_layout(rect=(0, 0, 1, .90))
    p = fig_dir / "4_machbarkeitskette.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    return p

def fig_energiebilanz(ts, rows, cyc, fig_dir, plt):
    yr = ts["days"] / 365.25
    fig, ax = plt.subplots(figsize=(11.5, 4.4))
    for t0, t1, _ in cyc["charge_intervals"]:
        ax.axvspan(t0 / YEAR, t1 / YEAR, color=C_ORANGE, alpha=.06, lw=0)
    ax.plot(yr, ts["E_in"], lw=1.6, color=C_ORANGE, label="eingespeichert (kumul.)")
    ax.plot(yr, ts["E_out"], lw=1.6, color=C_BLUE, label="entnommen (kumul.)")
    ax.plot(yr, ts["E_aq"], lw=1.2, color=C_AQUA, label="im Aquifer")
    ax.plot(yr, ts["E_cr"], lw=1.2, color=C_YELLOW, label="ins Deckgestein")
    if cyc["monthly_power_W"] and rows:
        dem = rows[0]["E_bedarf_GJ"]
        if dem:
            ax.plot(yr, dem * yr, lw=1.1, color=C_INK3, ls="--",
                    label="Bedarf (kumul.)")
    tot = ts["E_aq"][-1] + ts["E_cr"][-1]
    if tot > 0:
        ax.text(.012, .96, f"{100*ts['E_cr'][-1]/tot:.0f} % der gespeicherten "
                "Wärme liegt am Laufende im Deckgestein",
                transform=ax.transAxes, fontsize=9.5, color=C_INK, va="top",
                bbox=dict(fc="white", ec=C_RULE, lw=.8, pad=4))
    _style(ax, "Wohin geht die Wärme?",
           "Der Abstand zwischen „entnommen\" und „Bedarf\" ist der Deckungsgrad "
           "als Fläche.", "Betriebsjahr", "Energie [GJ]")
    ax.legend(fontsize=8.5, frameon=False, loc="upper left",
              bbox_to_anchor=(0, .88))
    fig.tight_layout()
    p = fig_dir / "5_energiebilanz.png"
    fig.savefig(p, dpi=140); plt.close(fig)
    return p

def _sample_plane(mesh, xs, ys, const, plane):
    """Feld auf ein regelmaessiges Raster in einer Ebene abtasten."""
    import pyvista as pv
    A, B = np.meshgrid(xs, ys)
    C = np.full_like(A, float(const))
    if plane == "xy":
        grid = pv.StructuredGrid(A, B, C)
    elif plane == "xz":
        grid = pv.StructuredGrid(A, C, B)
    else:
        grid = pv.StructuredGrid(A, B, C)
    res = grid.sample(mesh)
    T = np.asarray(res["T"]).reshape(A.shape, order="F") - 273.15
    valid = np.asarray(res["vtkValidPointMask"]).reshape(A.shape, order="F").astype(bool)
    T[~valid] = np.nan
    return A, B, T

def _pick_times(ts, cyc, n=6, whole=False):
    """n Zeitpunkte aus dem LETZTEN vollen Betriebsjahr, gleichmaessig verteilt."""
    t = ts["t"]
    if len(t) < 2:
        return []
    if whole:
        a, b = 0.0, t[-1]
    else:
        yr_last = max(0.0, np.floor(t[-1] / YEAR) - 1) if t[-1] >= 2 * YEAR else 0.0
        a = yr_last * YEAR
        b = min(a + YEAR, t[-1])
    want = np.linspace(a, b, n + 1)[1:]
    return [int(np.argmin(np.abs(t - w))) for w in want]

def fig_feldschnitte(cfg, run, geo, ts, cyc, fig_dir, plt):
    """Temperaturfeld zu mehreren Zeitpunkten - fixe Farbskala."""
    import pyvista as pv
    flow = bool(cfg.get("regional_gw", {}).get("enable", False))
    idx = _pick_times(ts, cyc, whole=flow)
    if not idx:
        return []
    T_amb = cyc["T_amb_K"] - 273.15
    T_inj = cfg["operation"]["T_hot_K"] - 273.15
    z0, z1 = geo.z_aq
    z_mid = 0.5 * (z0 + z1)
    lv = np.linspace(T_amb, T_inj, 26)
    out = []

    def panel(ax, A, B, T, title):
        cf = ax.contourf(A, B, T, levels=lv, cmap="inferno", extend="both")
        ax.contour(A, B, T, levels=[T_amb + 1.0], colors="#2a78d6", linewidths=.8)
        ax.contour(A, B, T, levels=[T_amb + 40.0], colors="white", linewidths=.9)
        ax.set_title(title, fontsize=9.5, color=C_INK, loc="left")
        ax.tick_params(labelsize=8, colors=C_INK2)
        return cf

    if run.axisym:
        r_max = max(30.0, 1.25 * float(np.nanmax(ts["r_front"])))
        rs = np.linspace(0.0, r_max, 200)
        zs = np.linspace(z0 - 45, z1 + 45, 190)
        xs = np.concatenate([-rs[::-1], rs])
        fig, axes = plt.subplots(2, 3, figsize=(14.5, 6.4), sharex=True, sharey=True)
        for ax, i in zip(axes.ravel(), idx):
            m = pv.read(run.snapshots[i][1])
            _A, B_half, T_half = _sample_plane(m, rs, zs, 0.0, "rz")
            T = np.hstack([T_half[:, ::-1], T_half])
            A = np.tile(xs, (len(zs), 1))
            B = np.hstack([B_half[:, ::-1], B_half])
            cf = panel(ax, A, B, T, f"Tag {ts['days'][i]:.0f}")
            for z in (z0, z1):
                ax.axhline(z, color="white", lw=.9, alpha=.75)
            ax.axvspan(-0.6, 0.6, color="#2a78d6", alpha=.25, lw=0)
        for ax in axes[-1]:
            ax.set_xlabel("r [m]", fontsize=9, color=C_INK2)
        for ax in axes[:, 0]:
            ax.set_ylabel("z [m]", fontsize=9, color=C_INK2)
        fig.colorbar(cf, ax=axes, shrink=.85, label="T [°C]")
        fig.suptitle("Temperaturfeld im letzten Betriebsjahr — Schnitt (r, z)",
                     x=.008, ha="left", fontsize=13, fontweight="bold", color=C_INK)
        fig.text(.008, .945, "Feste Farbskala T_amb…T_inj über alle Bilder. Blau = "
                 "1-K-Front, weiß = 50 °C. Blauer Balken = Filterstrecke.",
                 fontsize=9, color=C_INK3, va="top")
        p = fig_dir / "6_feldschnitt.png"
        fig.savefig(p, dpi=140, bbox_inches="tight"); plt.close(fig)
        out.append(p)
        return out

    xw = float(np.average(geo.r_cell[geo.mask_well] * 0 + 0))
    dom = pv.read(run.out_dir / f"{run.prefix}_domain.vtu")
    cc = dom.cell_centers().points
    mid = np.asarray(dom.cell_data["MaterialIDs"])
    xw = float(cc[mid == 3, 0].mean()) if (mid == 3).any() else 0.0
    yw = float(cc[mid == 3, 1].mean()) if (mid == 3).any() else 0.0
    dn = float(np.nanmax(ts["x_front_down"])) if "x_front_down" in ts else 0.0
    up = float(np.nanmin(ts["x_front_up"])) if "x_front_up" in ts else 0.0
    x0 = xw + min(up, -40) * 1.3
    x1 = xw + max(dn, 40) * 1.25
    yr = max(60.0, 0.7 * (x1 - x0) * 0.35)
    xs = np.linspace(x0, x1, 380)
    ys = np.linspace(yw - yr, yw + yr, 170)
    zs = np.linspace(z0 - 35, z1 + 35, 150)

    for tag, (bb, plane, const, ylab, ttl) in {
        "7_draufsicht": ((ys, "xy", z_mid), None, None, "y [m]",
                         "Draufsicht auf halber Aquiferhöhe — wandert die Fahne?"),
        "8_laengsschnitt": ((zs, "xz", yw), None, None, "z [m]",
                            "Längsschnitt durch die Brunnenachse — Auftrieb und Drift"),
    }.items():
        bs, plane, const = bb
        fig, axes = plt.subplots(3, 2, figsize=(14.5, 8.2), sharex=True, sharey=True)
        for ax, i in zip(axes.ravel(), idx):
            m = pv.read(run.snapshots[i][1])
            A, B, T = _sample_plane(m, xs, bs, const, plane)
            cf = panel(ax, A, B, T, f"Tag {ts['days'][i]:.0f}")
            ax.plot([xw], [yw if plane == "xy" else z_mid], "o", ms=5,
                    mfc="none", mec="#2a78d6", mew=1.6)
            if plane == "xz":
                for z in (z0, z1):
                    ax.axhline(z, color="white", lw=.9, alpha=.75)
        for ax in axes[-1]:
            ax.set_xlabel("x [m]  (Strömung nach +x)", fontsize=9, color=C_INK2)
        for ax in axes[:, 0]:
            ax.set_ylabel(ylab, fontsize=9, color=C_INK2)
        fig.colorbar(cf, ax=axes, shrink=.85, label="T [°C]")
        fig.suptitle(ttl, x=.008, ha="left", fontsize=13, fontweight="bold",
                     color=C_INK)
        fig.text(.008, .952, "Feste Farbskala T_amb…T_inj. Blau = 1-K-Front, "
                 "weiß = 50 °C. Kreis = Brunnen.",
                 fontsize=9, color=C_INK3, va="top")
        p = fig_dir / f"{tag}.png"
        fig.savefig(p, dpi=140, bbox_inches="tight"); plt.close(fig)
        out.append(p)
    return out

def schreibe_csv(cfg, run, geo, ts, cyc, mon, rows, chks, ziel):
    """Alle Zahlen des Laufs in EINEN Ordner, mit sprechenden Namen."""
    ziel = Path(ziel)
    ziel.mkdir(parents=True, exist_ok=True)
    K = 273.15
    geschrieben = []

    def _tab(name, kopf, zeilen):
        if not zeilen:
            return
        with open(ziel / name, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(kopf)
            w.writerows(zeilen)
        geschrieben.append(name)

    def _sp(name, key, f=None):
        a = ts.get(key)
        if a is None or not len(a):
            return None
        return (name, [(v if f is None else f(v)) for v in a])

    spalten = [x for x in (
        _sp("zeit_s",               "t"),
        _sp("tag",                  "days"),
        _sp("jahr",                 "days", lambda v: v / 365.25),
        _sp("T_brunnen_C",          "T_well", lambda v: v - K),
        _sp("T_sonde_aquifermitte_C", "T_probe", lambda v: v - K),
        _sp("T_max_C",              "T_max", lambda v: v - K),
        _sp("T_min_C",              "T_min", lambda v: v - K),
        _sp("T_injektion_C",        "T_inj", lambda v: v - K),
        _sp("massenstrom_kg_s",     "mdot"),
        _sp("brunnenleistung_kW",   "P", lambda v: v / 1e3),
        _sp("brunnendruck_kPa",     "dp_well", lambda v: v / 1e3),
        _sp("E_eingespeichert_GJ",  "E_in"),
        _sp("E_gefoerdert_GJ",      "E_out"),
        _sp("E_im_aquifer_GJ",      "E_aq"),
        _sp("E_im_deckgestein_GJ",  "E_cr"),
        _sp("fahne_radius_m",       "r_front"),
        _sp("fahne_stromab_m",      "x_front_down"),
        _sp("fahne_stromauf_m",     "x_front_up"),
        _sp("v_darcy_max_m_d",      "v_max", lambda v: v * 86400.0),
        _sp("aquifer_ueber_50C_pct", "frac_hot"),
    ) if x is not None]
    if spalten:
        n = len(spalten[0][1])
        zeilen = []
        for i in range(n):
            zeile = []
            for _, werte in spalten:
                v = werte[i]
                zeile.append("" if v != v else f"{v:.6g}")
            zeilen.append(zeile)
        _tab("zeitreihe.csv", [k for k, _ in spalten], zeilen)

    _tab("monatsbilanz.csv",
         ["monat", "betriebsart", "P_soll_kW", "T_brunnen_C",
          "E_gefordert_GJ", "E_geliefert_GJ", "deckung_pct"],
         [[m["monat"], "laden" if m["laden"] else "foerdern",
           f"{m['P_soll_kW']:.3f}", f"{m['T_well_C']:.3f}",
           f"{m['E_gefordert_GJ']:.4f}", f"{m['E_geliefert_GJ']:.4f}",
           f"{m['deckung_pct']:.2f}"] for m in (mon or [])])

    if rows:
        with open(ziel / "kennzahlen_jahr.csv", "w", newline="",
                  encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows([{k: (round(v, 4) if isinstance(v, float) else v)
                          for k, v in r.items()} for r in rows])
        geschrieben.append("kennzahlen_jahr.csv")

    _tab("pruefblatt.csv",
         ["pruefgroesse", "wert", "einheit", "ok_von", "ok_bis",
          "status", "diagnose"],
         [[c["name"], f"{c['val']:.4g}", c["unit"], c["band"][0],
           c["band"][1], c["status"], c["diag"]] for c in chks])

    def _flach(d, pre=""):
        o = []
        for k, v in d.items():
            if isinstance(v, dict):
                o += _flach(v, pre + k + ".")
            else:
                o.append((pre + k, v))
        return o
    _tab("konfiguration.csv", ["schluessel", "wert"],
         [[k, v] for k, v in _flach(cfg)])

    print(f"  [csv] {len(geschrieben)} Dateien -> {ziel}")
    return geschrieben

def report(cfg, out_dir=None, curves=None, rate_mult=None, report_dir=None):
    try:
        import sys as _s
        _s.stdout.reconfigure(errors="replace")
    except Exception:
        pass
    try:
        import pyvista as pv           # noqa: F401
    except ImportError:
        print("  [report] pyvista fehlt - uebersprungen.")
        return {"status": "no_pyvista"}

    run = detect_run(cfg, out_dir)
    if not run.snapshots:
        print(f"  [report] keine Ergebnisse in {run.out_dir}/ - uebersprungen.")
        return {"status": "no_results"}

    cyc = cycle_info(cfg, curves, rate_mult)
    _dom, geo = load_geometry(cfg, run)
    ts = build_timeseries(cfg, run, geo, cyc)
    rows, e_dem = energy_metrics(ts, cyc)
    mon = monthly_balance(ts, cyc, (rows[-1]["jahr"] - 1) if rows else 0)
    chks = checks(cfg, run, geo, ts, cyc, rows)

    rep = Path(report_dir) if report_dir else run.out_dir
    rep.mkdir(parents=True, exist_ok=True)
    fig_dir = rep / "figures"
    fig_dir.mkdir(exist_ok=True)

    schreibe_csv(cfg, run, geo, ts, cyc, mon, rows, chks, rep / "csv")

    print("\n  " + "=" * 74)
    print("  PRUEFBLATT" + (f"   ({run.prefix}, {ts['days'][-1]/365.25:.2f} Jahre "
                            f"gerechnet)" if len(ts["days"]) else ""))
    print("  " + "-" * 74)
    for c in chks:
        flag = {"OK": "  ok  ", "WARNUNG": " WARN ", "FEHLER": " FEHLER"}[c["status"]]
        print(f"  [{flag}] {c['name']:<34s} {c['val']:>11.4g} {c['unit']:<5s}"
              f" (ok: {c['band'][0]:g} .. {c['band'][1]:g})")
        if c["status"] != "OK":
            print(f"           -> {c['diag']}")
    if rows:
        print("  " + "-" * 74)
        hdr = f"  {'Jahr':>4} {'eta[%]':>7} {'etaF[%]':>8} "
        hdr += f"{'Deckung[%]':>11} " if e_dem else ""
        print(hdr + f"{'T_foerd[C]':>11} {'Deckgest.[GJ]':>14}")
        for r in rows:
            line = f"  {r['jahr']:>4} {r['eta_pct']:>7.1f} {r['eta_feld_pct']:>8.1f} "
            line += f"{r['deckung_pct']:>11.1f} " if e_dem else ""
            print(line + f"{r['T_foerder_C']:>11.2f} {r['E_deckgestein_GJ']:>14.0f}")
    print("  " + "=" * 74)

    figs = []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams.update({"font.size": 10, "figure.facecolor": "white",
                             "axes.facecolor": "white"})
    except ImportError:
        print("  [report] matplotlib fehlt - nur CSV geschrieben.")
        return {"status": "ok", "rows": rows, "checks": chks, "figures": []}

    for fn, args in (
        (fig_pruefblatt, (chks, fig_dir, plt)),
        (fig_brunnentemperatur, (cfg, ts, cyc, fig_dir, plt)),
        (fig_deckungsgrad, (rows, cyc, fig_dir, plt)),
        (fig_monatsbilanz, (mon, cyc, fig_dir, plt)),
        (fig_machbarkeit, (cfg, ts, cyc, mon, fig_dir, plt)),
        (fig_energiebilanz, (ts, rows, cyc, fig_dir, plt)),
        (fig_feldschnitte, (cfg, run, geo, ts, cyc, fig_dir, plt)),
    ):
        try:
            p = fn(*args)
            if isinstance(p, (list, tuple)):
                figs.extend([q for q in p if q])
            elif p:
                figs.append(p)
        except Exception as e:
            print(f"  [report] {fn.__name__}: {type(e).__name__}: {e}")
    print(f"  [report] {len(figs)} Abbildungen -> {fig_dir}")
    return {"status": "ok", "rows": rows, "checks": chks, "figures": figs,
            "timeseries": ts}

def auto_report(cfg, out_dir=None, curves=None, rate_mult=None, report_dir=None):
    """Wie report(), aber schluckt jeden Fehler: kein Lauf scheitert am Bericht."""
    if os.environ.get("ATES_REPORT", "1") == "0":
        return None
    try:
        return report(cfg, out_dir, curves, rate_mult, report_dir)
    except Exception as e:
        print(f"  [report] uebersprungen ({type(e).__name__}: {e})")
        return None

import sys as _sys
import types as _types
_mod = _types.ModuleType("ates_report")
for _n in ("auto_report", "report", "schreibe_csv"):
    if _n in globals():
        setattr(_mod, _n, globals()[_n])
_sys.modules["ates_report"] = _mod
#!/usr/bin/env python3
"""
ATES 3D demo for OpenGeoSys 6 (Hydro-Thermal process).

Geschichtetes Reservoir:  Cap Rock (oben) | Aquifer | Cap Rock (unten)
Ein Brunnenfilter als kleine Box im Aquifer (Single-Well-Anlage).
Zyklischer Lade-/Entlade-Betrieb über CurveScaled-Parameter:
  - Pressure-Equation: Volumetric Source ±Q/V_well auf Brunnenbox
  - Temperature: Dirichlet-BC auf Brunnenbox (curve-skaliert)
Der Lateralrand des Aquifers ist ein p=0-Outlet, damit das injizierte
Wasser entweichen kann.

Alle Modellgrößen sind im CONFIG-Block einstellbar.
Aufruf:
    python ates_3d.py            # Mesh + .prj + OGS-Lauf
    python ates_3d.py --no-run   # nur Setup, kein OGS
    python ates_3d.py --no-mesh  # nur .prj (Mesh muss existieren)
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

import gmsh
import re
import numpy as np

CONFIG: dict = {
    "domain": {
        "size_x_m":   400.0,
        "size_y_m":   250.0,
        "z_base_m":     0.0,
    },
    "layers": {
        "caprock_bottom_thickness_m": 40.0,
        "aquifer_thickness_m":        20.0,
        "caprock_top_thickness_m":    40.0,
    },
    "wells": {
        "hot_well_xy":   ( 0.0,  0.0),
        "screen_bottom_offset_m": 0.0,
        "screen_top_offset_m":    0.0,
        "screen_dx_m":             1.0,
        "screen_dy_m":             1.0,
        "screen_permeability_m2":  1.0e-11,
        "production_control":      "fixed",
        "max_rate_factor":         6.0,
        "demand_iterations":       3,
    },
    "regional_gw": {
        "enable":            True,
        "gradient_m_per_m":  2.0e-3,
        "direction_deg":     0.0,
    },
    "mesh": {
        "size_in_well_m":       0.5,
        "size_near_wells_m":    1.7,
        "size_far_m":          20.0,
        "well_size_radius_m":  22.0,
        "well_size_radius_far_m": 65.0,
    },
    "materials": {
        "aquifer": {
            "permeability_m2": 1.0e-11,
            "porosity":        0.25,
            "rho_s_kg_m3":  2650.0,
            "cp_s_J_kgK":    800.0,
            "lambda_s_W_mK":   3.0,
        },
        "caprock_top": {
            "permeability_m2": 1.0e-16,
            "porosity":        0.25,
            "rho_s_kg_m3":  2650.0,
            "cp_s_J_kgK":    800.0,
            "lambda_s_W_mK":   3.0,
        },
        "caprock_bottom": {
            "permeability_m2": 1.0e-16,
            "porosity":        0.25,
            "rho_s_kg_m3":  2650.0,
            "cp_s_J_kgK":    800.0,
            "lambda_s_W_mK":   3.0,
        },
    },
    "fluid": {
        "rho_ref_kg_m3":  1000.0,
        "T_ref_K":         283.15,
        "beta_1_per_K":   -4.0e-4,
        "viscosity_Pa_s":       1.3e-3,
        "visc_slope_1_per_K":  -1.28e-2,
        "cp_J_kgK":        4180.0,
        "lambda_W_mK":     0.6,
    },
    "dispersion": {
        "alpha_L_m": 1.0,
        "alpha_T_m": 0.1,
    },
    "initial": {
        "T_K":  283.15,
        "p_Pa": 1.0e5,
    },
    "operation": {
        "mass_flow_rate_kg_s": 3.0,
        "T_hot_K":  333.15,
        "T_cold_K": 283.15,
        "fluid_storage_1_per_Pa": 4.5e-10,
        "solid_storage_1_per_Pa": 1.0e-10,
    },
    "cycles": {
        "n_cycles":                     2,
        "charge_days":                 90,
        "storage_after_charge_days":    0,
        "discharge_days":              90,
        "storage_after_discharge_days": 0,
        "ramp_days":                    3.0,
        "monthly_power_W": [
            -400_000, -350_000, -120_000, 0,
            +250_000, +450_000, +500_000, +400_000, +120_000,
            -180_000, -320_000, -350_000,
        ],
        "monthly_T_inj_K":              None,
    },
    "time": {
        "dt_seconds":           86400.0,
        "output_every_n_steps": 5,
        "gravity":              True,
    },
    "output": {
        "prefix":    "ates_3d",
        "out_dir":   "out",
        "variables": ["T", "p", "darcy_velocity"],
    },
    "solver": {
        "solver_type":  "BiCGSTAB",
        "precon_type":  "ILUT",
        "linear_tol":   1.0e-10,
        "linear_iter":  10000,
        "nonlinear_iter": 20,
        "rel_tol_T":   1.0e-4,
        "rel_tol_p":   1.0e-4,
    },
}

DAY = 86400.0
G   = 9.81

def build_mesh(cfg: dict, out_dir: Path) -> Path:
    """Schichtmodell mit gmsh: 3 Schichten + 1 Brunnenbox im Aquifer."""
    msh_path = out_dir / f"{cfg['output']['prefix']}.msh"

    Lx     = cfg["domain"]["size_x_m"]
    Ly     = cfg["domain"]["size_y_m"]
    z_base = cfg["domain"]["z_base_m"]

    t_cb = cfg["layers"]["caprock_bottom_thickness_m"]
    t_aq = cfg["layers"]["aquifer_thickness_m"]
    t_ct = cfg["layers"]["caprock_top_thickness_m"]

    z_aq_bot = z_base + t_cb
    z_aq_top = z_aq_bot + t_aq
    z_top    = z_aq_top + t_ct

    hw  = cfg["wells"]["hot_well_xy"]
    sob = cfg["wells"]["screen_bottom_offset_m"]
    sot = cfg["wells"]["screen_top_offset_m"]
    dx  = cfg["wells"]["screen_dx_m"]
    dy  = cfg["wells"]["screen_dy_m"]
    z_screen_bot = z_aq_bot + sob
    z_screen_top = z_aq_top - sot
    h_screen = z_screen_top - z_screen_bot

    s_in   = cfg["mesh"]["size_in_well_m"]
    s_near = cfg["mesh"]["size_near_wells_m"]
    s_far  = cfg["mesh"]["size_far_m"]
    r_near = cfg["mesh"]["well_size_radius_m"]
    r_far  = cfg["mesh"]["well_size_radius_far_m"]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("ates")

    x0, y0 = -Lx / 2.0, -Ly / 2.0
    box_cb = gmsh.model.occ.addBox(x0, y0, z_base,   Lx, Ly, t_cb)
    box_aq = gmsh.model.occ.addBox(x0, y0, z_aq_bot, Lx, Ly, t_aq)
    box_ct = gmsh.model.occ.addBox(x0, y0, z_aq_top, Lx, Ly, t_ct)
    box_hw = gmsh.model.occ.addBox(hw[0] - dx / 2.0, hw[1] - dy / 2.0,
                                   z_screen_bot, dx, dy, h_screen)

    gmsh.model.occ.fragment([(3, box_cb)],
                            [(3, box_aq), (3, box_ct), (3, box_hw)])
    gmsh.model.occ.synchronize()

    vol_aq, vol_ct, vol_cb, vol_hw = [], [], [], []
    for dim, tag in gmsh.model.getEntities(3):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(dim, tag)
        zc = 0.5 * (zmin + zmax)
        ext_x = xmax - xmin
        if zc > z_aq_top + 1e-6:
            vol_ct.append(tag)
        elif zc < z_aq_bot - 1e-6:
            vol_cb.append(tag)
        else:
            xc = 0.5 * (xmin + xmax)
            yc = 0.5 * (ymin + ymax)
            small = ext_x < 0.5 * Lx
            if small and abs(xc - hw[0]) < dx and abs(yc - hw[1]) < dy:
                vol_hw.append(tag)
            else:
                vol_aq.append(tag)
    if not vol_aq or not vol_hw:
        raise RuntimeError("Volumenklassifizierung fehlgeschlagen (gmsh-Fragmentierung).")

    surf_top, surf_bot, surf_lat_aq = [], [], []
    for dim, tag in gmsh.model.getEntities(2):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(dim, tag)
        zc = 0.5 * (zmin + zmax)
        if (xmax - xmin) >= 0.9 * Lx and abs(zc - z_top) < 1e-6:
            surf_top.append(tag)
            continue
        if (xmax - xmin) >= 0.9 * Lx and abs(zc - z_base) < 1e-6:
            surf_bot.append(tag)
            continue
        on_x_edge = abs(xmin - x0) < 1e-6 and abs(xmax - x0) < 1e-6
        on_x_edge_pos = abs(xmin - (x0 + Lx)) < 1e-6 and abs(xmax - (x0 + Lx)) < 1e-6
        on_y_edge = abs(ymin - y0) < 1e-6 and abs(ymax - y0) < 1e-6
        on_y_edge_pos = abs(ymin - (y0 + Ly)) < 1e-6 and abs(ymax - (y0 + Ly)) < 1e-6
        on_outer = on_x_edge or on_x_edge_pos or on_y_edge or on_y_edge_pos
        in_aquifer_z = (zmin >= z_aq_bot - 1e-6) and (zmax <= z_aq_top + 1e-6)
        if on_outer and in_aquifer_z:
            surf_lat_aq.append(tag)

    surf_hw = []
    for tag in vol_hw:
        for d, t in gmsh.model.getBoundary([(3, tag)], oriented=False):
            if d == 2:
                surf_hw.append(abs(t))
    surf_hw = sorted(set(surf_hw))

    gmsh.model.addPhysicalGroup(3, vol_aq, tag=1, name="aquifer")
    gmsh.model.addPhysicalGroup(3, vol_ct, tag=2, name="caprock_top")
    gmsh.model.addPhysicalGroup(3, vol_cb, tag=3, name="caprock_bottom")
    gmsh.model.addPhysicalGroup(3, vol_hw, tag=4, name="hot_well_vol")
    gmsh.model.addPhysicalGroup(2, surf_top, tag=10, name="top")
    gmsh.model.addPhysicalGroup(2, surf_bot, tag=11, name="bottom")
    gmsh.model.addPhysicalGroup(2, surf_hw,  tag=12, name="hot_well_surf")
    if surf_lat_aq:
        gmsh.model.addPhysicalGroup(2, surf_lat_aq, tag=14, name="lateral_aquifer")

    _gw = cfg.get("regional_gw", {})
    if surf_lat_aq and _gw.get("enable", False):
        import math as _m
        _ang = _m.radians(float(_gw.get("direction_deg", 0.0)))
        _ex, _ey = _m.cos(_ang), _m.sin(_ang)
        _proj = {}
        for _t in surf_lat_aq:
            _bb = gmsh.model.getBoundingBox(2, _t)
            _cx, _cy = 0.5 * (_bb[0] + _bb[3]), 0.5 * (_bb[1] + _bb[4])
            _proj[_t] = _cx * _ex + _cy * _ey
        _lo, _hi = min(_proj.values()), max(_proj.values())
        _thr = _lo + 0.25 * (_hi - _lo)
        _inflow = sorted(t for t, v in _proj.items() if v <= _thr)
        if _inflow:
            gmsh.model.addPhysicalGroup(2, _inflow, tag=15,
                                        name="lateral_inflow")

    well_surfaces: list[int] = []
    for tag in vol_hw:
        for d, t in gmsh.model.getBoundary([(3, tag)], oriented=False):
            if d == 2:
                well_surfaces.append(abs(t))

    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "SurfacesList", list(set(well_surfaces)))
    f_thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_thr, "InField", f_dist)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMin", s_near)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMax", s_far)
    gmsh.model.mesh.field.setNumber(f_thr, "DistMin", r_near)
    gmsh.model.mesh.field.setNumber(f_thr, "DistMax", r_far)
    gmsh.model.mesh.field.setAsBackgroundMesh(f_thr)

    well_points: list[tuple[int, int]] = []
    for tag in vol_hw:
        for d, t in gmsh.model.getBoundary([(3, tag)], recursive=True, oriented=False):
            if d == 0:
                well_points.append((d, t))
    if well_points:
        gmsh.model.mesh.setSize(well_points, s_in)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

    gmsh.model.mesh.generate(3)
    gmsh.write(str(msh_path))
    gmsh.finalize()
    return msh_path

def _mesh_files(cfg: dict) -> dict[str, str]:
    prefix = cfg["output"]["prefix"]
    d = {
        "domain":          f"{prefix}_domain.vtu",
        "top":             f"{prefix}_physical_group_top.vtu",
        "bottom":          f"{prefix}_physical_group_bottom.vtu",
        "hot_well_vol":    f"{prefix}_physical_group_hot_well_vol.vtu",
        "hot_well_surf":   f"{prefix}_physical_group_hot_well_surf.vtu",
        "lateral_aquifer": f"{prefix}_physical_group_lateral_aquifer.vtu",
    }
    if cfg.get("regional_gw", {}).get("enable", False):
        d["lateral_inflow"] = f"{prefix}_physical_group_lateral_inflow.vtu"
    return d

def _safe_name(name):
    """Physical-Group-Namen plattform- und dateisystemsicher machen."""
    keep = ("abcdefghijklmnopqrstuvwxyz"
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")
    return "".join(c if c in keep else "_" for c in str(name))

def convert_mesh(cfg: dict, msh_path: Path, out_dir: Path) -> dict[str, str]:
    """gmsh-Mesh -> OGS .vtu (Domäne + Subdomänen)."""
    import ogstools as ot
    prefix = cfg["output"]["prefix"]
    meshes = ot.Meshes.from_gmsh(
        filename=str(msh_path), dim=3, reindex=True, log=False
    )
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, mesh in meshes.items():
        if name == "domain":
            fname = f"{prefix}_domain.vtu"
        else:
            fname = f"{prefix}_physical_group_{_safe_name(name)}.vtu"
        mesh.save(str(out_dir / fname), binary=True)
    return _mesh_files(cfg)

def build_cycle_curves(cfg: dict, rate_mult=None) -> dict:
    """Normierte Pump-/Leistungskurve g(t) ∈ [−1,1] + Beladungs-Intervalle"""
    cyc  = cfg["cycles"]
    n    = cyc["n_cycles"]
    ramp = max(60.0, cyc["ramp_days"] * DAY)
    T0     = cfg["initial"]["T_K"]
    T_hot  = cfg["operation"]["T_hot_K"]
    cp_f   = cfg["fluid"]["cp_J_kgK"]
    dT_ref = T_hot - T0
    if dT_ref <= 0:
        raise ValueError("T_hot_K muss > initial.T_K sein (ΔT_ref > 0).")

    def _piecewise(seg_values, seg_durs):
        times = [0.0]; vals = [0.0]; t_now = 0.0
        for _ in range(n):
            for g, dur in zip(seg_values, seg_durs):
                if dur <= 0:
                    continue
                t_now += ramp
                times.append(t_now); vals.append(g)
                hold = max(0.0, dur - ramp)
                if hold > 0:
                    t_now += hold
                    times.append(t_now); vals.append(g)
        t_now += ramp
        times.append(t_now); vals.append(0.0)
        return np.array(times), np.array(vals), t_now

    def _merge(ivs):
        m = []
        for t0, t1, Ti in ivs:
            if m and abs(m[-1][1] - t0) < 1e-6 and m[-1][2] == Ti:
                m[-1] = (m[-1][0], t1, Ti)
            else:
                m.append((t0, t1, Ti))
        return m

    monthly_P = cyc.get("monthly_power_W")
    if monthly_P is not None:
        assert len(monthly_P) == 12, "cycles.monthly_power_W muss 12 Werte enthalten."
        monthly_T = cyc.get("monthly_T_inj_K")
        P = np.asarray(monthly_P, dtype=float)
        P_nom = float(np.max(np.abs(P)))
        if P_nom == 0.0:
            raise ValueError("cycles.monthly_power_W enthält nur Nullen.")
        mdot_nom  = P_nom / (cp_f * dT_ref)
        month_dur = 365.25 / 12.0 * DAY
        g_vals = (P / P_nom).tolist()
        if rate_mult is not None:
            g_vals = [g * (rate_mult[m] if P[m] < 0 else 1.0) for m, g in enumerate(g_vals)]
        times, vals, t_tot = _piecewise(g_vals, [month_dur] * 12)
        ivs = []
        for y in range(n):
            for m in range(12):
                if P[m] > 0.0:
                    k = y * 12 + m
                    Ti = float(monthly_T[m]) if monthly_T is not None else T_hot
                    ivs.append((k * month_dur, (k + 1) * month_dur, Ti))
        return {"t_total": t_tot, "cycle_power": (times, vals),
                "charge_intervals": _merge(ivs), "P_nom_W": P_nom,
                "mdot_nom_kg_s": mdot_nom, "mode": "monthly"}

    mdot_nom = cfg["operation"]["mass_flow_rate_kg_s"]
    P_nom    = mdot_nom * cp_f * dT_ref
    durs = [cyc["charge_days"] * DAY, cyc["storage_after_charge_days"] * DAY,
            cyc["discharge_days"] * DAY, cyc["storage_after_discharge_days"] * DAY]
    times, vals, t_tot = _piecewise([+1.0, 0.0, -1.0, 0.0], durs)
    cyc_dur = sum(durs)
    charge_intervals = [(y * cyc_dur, y * cyc_dur + durs[0], T_hot)
                        for y in range(n) if durs[0] > 0]
    return {"t_total": t_tot, "cycle_power": (times, vals),
            "charge_intervals": charge_intervals, "P_nom_W": P_nom,
            "mdot_nom_kg_s": mdot_nom, "mode": "phases"}

def _se(parent: ET.Element, tag: str, text=None, **attrs) -> ET.Element:
    el = ET.SubElement(parent, tag, **{k: str(v) for k, v in attrs.items()})
    if text is not None:
        el.text = str(text)
    return el

def _add_const_property(parent: ET.Element, name: str, value: float) -> None:
    p = _se(parent, "property")
    _se(p, "name", name)
    _se(p, "type", "Constant")
    _se(p, "value", value)

def _add_phase_aqueous(phases: ET.Element, fluid: dict, op: dict) -> None:
    ph = _se(phases, "phase")
    _se(ph, "type", "AqueousLiquid")
    props = _se(ph, "properties")

    p = _se(props, "property")
    _se(p, "name", "density")
    beta = fluid.get("beta_1_per_K", 0.0)
    if abs(beta) > 1e-12:
        _se(p, "type", "Linear")
        _se(p, "reference_value", fluid["rho_ref_kg_m3"])
        iv = _se(p, "independent_variable")
        _se(iv, "variable_name", "temperature")
        _se(iv, "reference_condition", fluid["T_ref_K"])
        _se(iv, "slope", fluid["beta_1_per_K"])
    else:
        _se(p, "type", "Constant")
        _se(p, "value", fluid["rho_ref_kg_m3"])

    vslope = fluid.get("visc_slope_1_per_K")
    if vslope:
        pv = _se(props, "property"); _se(pv, "name", "viscosity"); _se(pv, "type", "Linear")
        _se(pv, "reference_value", fluid["viscosity_Pa_s"])
        iv = _se(pv, "independent_variable")
        _se(iv, "variable_name", "temperature")
        _se(iv, "reference_condition", fluid["T_ref_K"])
        _se(iv, "slope", vslope)
    else:
        _add_const_property(props, "viscosity", fluid["viscosity_Pa_s"])
    _add_const_property(props, "specific_heat_capacity", fluid["cp_J_kgK"])
    _add_const_property(props, "thermal_conductivity",   fluid["lambda_W_mK"])
    _add_const_property(props, "storage",                op["fluid_storage_1_per_Pa"])

def _add_phase_solid(phases: ET.Element, m: dict, op: dict) -> None:
    ph = _se(phases, "phase")
    _se(ph, "type", "Solid")
    props = _se(ph, "properties")
    _add_const_property(props, "density",                m["rho_s_kg_m3"])
    _add_const_property(props, "specific_heat_capacity", m["cp_s_J_kgK"])
    _add_const_property(props, "thermal_conductivity",   m["lambda_s_W_mK"])
    _add_const_property(props, "storage",                op["solid_storage_1_per_Pa"])

def _add_medium(media: ET.Element, mid: int, mat: dict, fluid: dict,
                op: dict, disp: dict) -> None:
    med = _se(media, "medium", id=mid)
    phases = _se(med, "phases")
    _add_phase_aqueous(phases, fluid, op)
    _add_phase_solid(phases, mat, op)
    props = _se(med, "properties")
    _add_const_property(props, "porosity",     mat["porosity"])
    _add_const_property(props, "permeability", mat["permeability_m2"])

    p = _se(props, "property")
    _se(p, "name", "thermal_conductivity")
    _se(p, "type", "EffectiveThermalConductivityPorosityMixing")

    _add_const_property(props, "thermal_longitudinal_dispersivity", disp["alpha_L_m"])
    _add_const_property(props, "thermal_transversal_dispersivity",  disp["alpha_T_m"])
    _add_const_property(props, "storage", 0.0)

def _curve_xml(parent: ET.Element, name: str, t: np.ndarray, v: np.ndarray) -> None:
    c = _se(parent, "curve")
    _se(c, "name", name)
    _se(c, "coords", " ".join(f"{x:.6e}" for x in t))
    _se(c, "values", " ".join(f"{x:.6e}" for x in v))

def _curve_scaled_param(parent: ET.Element, name: str,
                        curve_name: str, base_param_name: str) -> None:
    p = _se(parent, "parameter")
    _se(p, "name", name)
    _se(p, "type", "CurveScaled")
    _se(p, "curve", curve_name)
    _se(p, "parameter", base_param_name)

def _const_param(parent: ET.Element, name: str, value: float) -> None:
    p = _se(parent, "parameter")
    _se(p, "name", name)
    _se(p, "type", "Constant")
    _se(p, "value", value)

def build_prj(cfg: dict, out_dir: Path, mesh_files: dict[str, str], curves: dict) -> Path:
    prefix = cfg["output"]["prefix"]
    fluid  = cfg["fluid"]
    op     = cfg["operation"]
    init   = cfg["initial"]
    sol    = cfg["solver"]

    h_screen = (cfg["layers"]["aquifer_thickness_m"]
                - cfg["wells"]["screen_top_offset_m"]
                - cfg["wells"]["screen_bottom_offset_m"])
    dx_w = cfg["wells"]["screen_dx_m"]
    dy_w = cfg["wells"]["screen_dy_m"]
    V_well = dx_w * dy_w * h_screen

    rho_f = fluid["rho_ref_kg_m3"]
    q_mass_amp = curves["mdot_nom_kg_s"] / rho_f / V_well

    root = ET.Element("OpenGeoSysProject")

    meshes = ET.SubElement(root, "meshes")
    _mesh_keys = ["domain", "top", "bottom", "hot_well_vol", "hot_well_surf"]
    if "lateral_aquifer" in mesh_files:
        _mesh_keys.append("lateral_aquifer")
    if "lateral_inflow" in mesh_files:
        _mesh_keys.append("lateral_inflow")
    for k in _mesh_keys:
        _se(meshes, "mesh", mesh_files[k])

    processes = _se(root, "processes")
    proc = _se(processes, "process")
    _se(proc, "name", "HT")
    _se(proc, "type", "HT")
    _se(proc, "integration_order", 2)
    pv = _se(proc, "process_variables")
    _se(pv, "temperature", "T")
    _se(pv, "pressure",    "p")
    sv = _se(proc, "secondary_variables")
    _se(sv, "secondary_variable", internal_name="darcy_velocity",
        output_name="darcy_velocity")
    bf = "0 0 -9.81" if cfg["time"]["gravity"] else "0 0 0"
    _se(proc, "specific_body_force", bf)

    well_mat = dict(cfg["materials"]["aquifer"])
    well_mat["permeability_m2"] = cfg["wells"]["screen_permeability_m2"]
    disp = cfg["dispersion"]
    media = _se(root, "media")
    _add_medium(media, 0, cfg["materials"]["aquifer"],        fluid, op, disp)
    _add_medium(media, 1, cfg["materials"]["caprock_top"],    fluid, op, disp)
    _add_medium(media, 2, cfg["materials"]["caprock_bottom"], fluid, op, disp)
    _add_medium(media, 3, well_mat,                           fluid, op, disp)

    tl = _se(root, "time_loop")
    procs = _se(tl, "processes")
    p_ref = _se(procs, "process", ref="HT")
    _se(p_ref, "nonlinear_solver", "basic_picard")
    cc = _se(p_ref, "convergence_criterion")
    _se(cc, "type",      "PerComponentDeltaX")
    _se(cc, "norm_type", "NORM2")
    _se(cc, "reltols",   f"{sol['rel_tol_T']} {sol['rel_tol_p']}")
    td = _se(p_ref, "time_discretization")
    _se(td, "type", "BackwardEuler")

    dt0 = cfg["time"]["dt_seconds"]
    ts = _se(p_ref, "time_stepping")
    _se(ts, "type",        "IterationNumberBasedTimeStepping")
    _se(ts, "t_initial",   0.0)
    _se(ts, "t_end",       curves["t_total"])
    _se(ts, "initial_dt",  dt0)
    _se(ts, "minimum_dt",  dt0 / 64.0)
    _se(ts, "maximum_dt",  dt0)
    _se(ts, "number_iterations", "1 4 8 12")
    _se(ts, "multiplier",        "1.5 1.0 0.5 0.25")

    out_step = dt0 * cfg["time"]["output_every_n_steps"]
    n_out    = int(np.ceil(curves["t_total"] / out_step))
    out_times = sorted(set(
        [0.0] + [min((k + 1) * out_step, curves["t_total"]) for k in range(n_out)]
    ))
    out = _se(tl, "output")
    _se(out, "type",   "VTK")
    _se(out, "prefix", prefix)
    _se(out, "fixed_output_times", " ".join(f"{t:.6f}" for t in out_times))
    _se(out, "output_iteration_results", "false")
    vars_el = _se(out, "variables")
    for v in cfg["output"]["variables"]:
        _se(vars_el, "variable", v)

    params = _se(root, "parameters")
    _const_param(params, "T0", init["T_K"])

    z_top = (cfg["domain"]["z_base_m"] + cfg["layers"]["caprock_bottom_thickness_m"]
             + cfg["layers"]["aquifer_thickness_m"] + cfg["layers"]["caprock_top_thickness_m"])
    rho0  = fluid["rho_ref_kg_m3"]; p_top = init["p_Pa"]
    hydro = (f"{p_top:.6g} + {rho0*9.81:.6g}*({z_top:.6g} - z)"
             if cfg["time"]["gravity"] else f"{p_top:.6g}")
    p_el = _se(params, "parameter"); _se(p_el, "name", "p0"); _se(p_el, "type", "Function")
    _se(p_el, "expression", hydro)
    gw = cfg.get("regional_gw", {})
    if gw.get("enable", False):
        import math
        rho_g = rho0 * 9.81
        i = float(gw["gradient_m_per_m"]); alpha = math.radians(float(gw["direction_deg"]))
        gx = -rho_g * i * math.cos(alpha); gy = -rho_g * i * math.sin(alpha)
        p_lat = _se(params, "parameter"); _se(p_lat, "name", "p_lateral"); _se(p_lat, "type", "Function")
        _se(p_lat, "expression", f"{hydro} + ({gx:.6g})*x + ({gy:.6g})*y")

    _const_param(params, "q_mass_amp", q_mass_amp)
    _curve_scaled_param(params, "q_mass_well", "cycle_power", "q_mass_amp")
    tinj_values = sorted({round(iv[2], 6) for iv in curves["charge_intervals"]})
    tinj_param = {}
    for idx, Ti in enumerate(tinj_values):
        nm = "T_inj" if len(tinj_values) == 1 else f"T_inj_{idx}"
        _const_param(params, nm, Ti); tinj_param[round(Ti, 6)] = nm

    cv = _se(root, "curves")
    t, v = curves["cycle_power"]; _curve_xml(cv, "cycle_power", t, v)

    pvars = _se(root, "process_variables")

    pv_T = _se(pvars, "process_variable")
    _se(pv_T, "name", "T")
    _se(pv_T, "components", 1)
    _se(pv_T, "order", 1)
    _se(pv_T, "initial_condition", "T0")
    bcs_T = _se(pv_T, "boundary_conditions")
    for face in ("top", "bottom"):
        bc = _se(bcs_T, "boundary_condition")
        _se(bc, "mesh",      Path(mesh_files[face]).stem)
        _se(bc, "type",      "Dirichlet")
        _se(bc, "parameter", "T0")
    if "lateral_inflow" in mesh_files:
        bc = _se(bcs_T, "boundary_condition")
        _se(bc, "mesh",      Path(mesh_files["lateral_inflow"]).stem)
        _se(bc, "type",      "Dirichlet")
        _se(bc, "parameter", "T0")
    well_mesh_T = Path(mesh_files["hot_well_vol"]).stem
    for t0, t1, Ti in curves["charge_intervals"]:
        bc = _se(bcs_T, "boundary_condition")
        _se(bc, "mesh", well_mesh_T)
        _se(bc, "type", "DirichletWithinTimeInterval")
        _se(bc, "parameter", tinj_param[round(Ti, 6)])
        ti = _se(bc, "time_interval")
        _se(ti, "start", f"{t0:.6e}"); _se(ti, "end", f"{t1:.6e}")

    pv_p = _se(pvars, "process_variable")
    _se(pv_p, "name", "p")
    _se(pv_p, "components", 1)
    _se(pv_p, "order", 1)
    _se(pv_p, "initial_condition", "p0")
    bcs_p = _se(pv_p, "boundary_conditions")
    for face in ("top", "bottom"):
        bc = _se(bcs_p, "boundary_condition")
        _se(bc, "mesh",      Path(mesh_files[face]).stem)
        _se(bc, "type",      "Dirichlet")
        _se(bc, "parameter", "p0")
    if "lateral_aquifer" in mesh_files:
        bc = _se(bcs_p, "boundary_condition")
        _se(bc, "mesh",      Path(mesh_files["lateral_aquifer"]).stem)
        _se(bc, "type",      "Dirichlet")
        if cfg.get("regional_gw", {}).get("enable", False):
            _se(bc, "parameter", "p_lateral")
        else:
            _se(bc, "parameter", "p0")
    sts_p = _se(pv_p, "source_terms")
    st = _se(sts_p, "source_term")
    _se(st, "mesh",      Path(mesh_files["hot_well_vol"]).stem)
    _se(st, "type",      "Volumetric")
    _se(st, "parameter", "q_mass_well")

    nls = _se(root, "nonlinear_solvers")
    n = _se(nls, "nonlinear_solver")
    _se(n, "name",          "basic_picard")
    _se(n, "type",          "Picard")
    _se(n, "max_iter",      sol["nonlinear_iter"])
    _se(n, "linear_solver", "general_linear_solver")

    lss = _se(root, "linear_solvers")
    ls = _se(lss, "linear_solver")
    _se(ls, "name", "general_linear_solver")
    eig = _se(ls, "eigen")
    stype = sol.get("solver_type", "BiCGSTAB")
    _se(eig, "solver_type", stype)
    if stype != "SparseLU":
        _se(eig, "precon_type",        sol.get("precon_type", "ILUT"))
        _se(eig, "max_iteration_step", sol["linear_iter"])
        _se(eig, "error_tolerance",    sol["linear_tol"])
    _se(eig, "scaling",            "true")

    _indent(root)
    prj_path = out_dir / f"{prefix}.prj"
    tree = ET.ElementTree(root)
    tree.write(prj_path, encoding="ISO-8859-1", xml_declaration=True)
    return prj_path

def _indent(elem: ET.Element, level: int = 0) -> None:
    pad = "\n" + level * "    "
    if len(elem):
        if not (elem.text and elem.text.strip()):
            elem.text = pad + "    "
        for child in elem:
            _indent(child, level + 1)
        if not (elem[-1].tail and elem[-1].tail.strip()):
            elem[-1].tail = pad
    if level and not (elem.tail and elem.tail.strip()):
        elem.tail = pad

def run_ogs(prj_path: Path) -> int:
    ogs_exe = shutil.which("ogs") or shutil.which("ogs.exe")
    if not ogs_exe:
        print("ogs.exe nicht im PATH gefunden – bitte installieren oder --no-run nutzen.",
              file=sys.stderr)
        return 1
    cmd = [ogs_exe, str(prj_path), "-o", str(prj_path.parent)]
    print(">>", " ".join(cmd))
    return subprocess.call(cmd)

def measure_produced_T(cfg: dict, out_dir: Path) -> dict:
    """|P|-gewichtete mittlere Brunnentemperatur je Fördermonat (P<0),"""
    import pyvista as pv
    prefix = cfg["output"]["prefix"]
    files = sorted(out_dir.glob(f"{prefix}_ts_*.vtu"),
                   key=lambda p: int(re.search(r"_ts_(\d+)_", p.name).group(1)))
    if not files:
        return {}
    z_mid = (cfg["domain"]["z_base_m"] + cfg["layers"]["caprock_bottom_thickness_m"]
             + cfg["layers"]["aquifer_thickness_m"] / 2.0)
    well_pt = pv.PolyData(np.array([[0.5, 0.0, z_mid]]))
    YEAR = 365.25; month_dur = YEAR / 12.0
    n_cyc = cfg["cycles"]["n_cycles"]; monthly_P = cfg["cycles"]["monthly_power_W"]
    y_start = max(1, n_cyc // 2)
    days, Tw = [], []
    for f in files:
        d = float(re.search(r"_t_([0-9.]+)\.vtu", f.name).group(1)) / DAY
        if d < y_start * YEAR:
            continue
        Tw.append(float(well_pt.sample(pv.read(f))["T"][0])); days.append(d)
    days = np.array(days); Tw = np.array(Tw); res = {}
    for mo in range(12):
        if monthly_P[mo] >= 0:
            continue
        vals = []
        for y in range(y_start, n_cyc):
            seg = (days >= y*YEAR + mo*month_dur) & (days < y*YEAR + (mo+1)*month_dur)
            if seg.any():
                vals.append(Tw[seg].mean())
        if vals:
            res[mo] = float(np.mean(vals))
    return res

def main() -> int:
    ap = argparse.ArgumentParser(description="ATES 3D OGS demo")
    ap.add_argument("--no-mesh", action="store_true", help="Mesh nicht neu erzeugen")
    ap.add_argument("--no-run",  action="store_true", help="OGS nicht ausführen")
    ap.add_argument("--years", type=int, default=None, help="Override cycles.n_cycles")
    ap.add_argument("--production-control", choices=["demand", "fixed"], default=None,
                    help="Override wells.production_control")
    args = ap.parse_args()
    if args.years is not None:
        CONFIG["cycles"]["n_cycles"] = args.years
    if args.production_control is not None:
        CONFIG["wells"]["production_control"] = args.production_control

    out_dir = Path(CONFIG["output"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = CONFIG["output"]["prefix"]
    msh_path = out_dir / f"{prefix}.msh"

    if not args.no_mesh:
        print("[1/3] gmsh: 3D Schichtmodell + Brunnenbox ...")
        build_mesh(CONFIG, out_dir)
        print(f"      {msh_path}")
        print("[2/3] msh2vtu: Konvertierung in OGS-Meshes ...")
        mesh_files = convert_mesh(CONFIG, msh_path, out_dir)
    else:
        mesh_files = _mesh_files(CONFIG)

    monthly = CONFIG["cycles"].get("monthly_power_W") is not None
    demand  = (CONFIG["wells"].get("production_control", "fixed") == "demand"
               and monthly and not args.no_run)
    T_amb  = CONFIG["initial"]["T_K"]
    dT_ref = CONFIG["operation"]["T_hot_K"] - T_amb
    max_fac = CONFIG["wells"].get("max_rate_factor", 6.0)
    rate_mult = [1.0] * 12
    n_iter = CONFIG["wells"].get("demand_iterations", 3) if demand else 1

    for it in range(n_iter):
        tag = f" (Bedarfs-Iteration {it+1}/{n_iter})" if demand else ""
        print(f"[3/3] OGS-Projektdatei erzeugen ...{tag}")
        curves = build_cycle_curves(CONFIG, rate_mult=rate_mult if demand else None)
        prj_path = build_prj(CONFIG, out_dir, mesh_files, curves)
        print(f"      {prj_path}  (t_end = {curves['t_total']/DAY:.1f} d)")
        if args.no_run:
            return 0
        print(">>> OGS starten")
        rc = run_ogs(prj_path)
        if rc != 0:
            return rc
        if demand and it < n_iter - 1:
            Tp = measure_produced_T(CONFIG, out_dir)
            for mo, T in Tp.items():
                target = min(dT_ref / max(T - T_amb, 1.0), max_fac)
                rate_mult[mo] = 0.5 * rate_mult[mo] + 0.5 * target
            print("  T_prod_avg [degC]:", {mo: round(T-273.15, 1) for mo, T in Tp.items()})
            print("  rate_factors:", [round(x, 2) for x in rate_mult])
    try:
        import sys as _sys
        _here = Path(__file__).resolve()
        for _p in (_here.parent, *_here.parents):
            if (_p / "ates_report.py").exists():
                _sys.path.insert(0, str(_p))
                break
        import ates_report
        ates_report.auto_report(CONFIG, out_dir, curves=curves,
                                rate_mult=(rate_mult if demand else None))
    except Exception as _e:
        print(f"  [report] uebersprungen: {type(_e).__name__}: {_e}")

    return 0

C = CONFIG
G = 9.81
RHO_F = C["fluid"]["rho_ref_kg_m3"]

T_AQ_K = FALL["T_aquifer_C"] + 273.15
C["initial"]["T_K"] = T_AQ_K
C["fluid"]["T_ref_K"] = T_AQ_K
C["operation"]["T_hot_K"] = FALL["T_injektion_C"] + 273.15

MU = C["fluid"]["viscosity_Pa_s"]

_T_NULL = FALL["T_aquifer_C"] - 1.0 / C["fluid"]["visc_slope_1_per_K"]
if FALL["T_injektion_C"] > _T_NULL - 10.0:
    raise ValueError(
        f"T_injektion_C = {FALL['T_injektion_C']:.0f} GradC ist zu hoch fuer "
        f"das lineare Viskositaetsmodell: mu(T) wird bei "
        f"{_T_NULL:.1f} GradC null und darueber negativ.\n"
        f"Zulaessig ist hier bis {_T_NULL - 10.0:.1f} GradC. Fuer hoehere "
        f"Injektionstemperaturen muss fluid.visc_slope_1_per_K angepasst "
        f"werden.")

def _k(kf_m_s: float) -> float:
    """Durchlaessigkeitsbeiwert kf [m/s] -> Permeabilitaet k [m2]."""
    return kf_m_s * MU / (RHO_F * G)

C["regional_gw"]["enable"] = True
C["regional_gw"]["gradient_m_per_m"] = FALL["gw_gradient"]
C["regional_gw"]["direction_deg"] = FALL["gw_richtung_grad"]

_a = FALL["aquifer"]
C["layers"]["aquifer_thickness_m"] = _a["maechtigkeit_m"]
aq = C["materials"]["aquifer"]
aq["permeability_m2"] = _k(_a["kf_m_s"])
aq["porosity"] = _a["porositaet"]
aq["rho_s_kg_m3"] = _a["dichte_korn_kg_m3"]
aq["cp_s_J_kgK"] = _a["waermekapazitaet_korn_J_kgK"]
aq["lambda_s_W_mK"] = _a["waermeleitfaehigkeit_korn_W_mK"]

_d = FALL["deckgestein"]
C["layers"]["caprock_top_thickness_m"] = _d["maechtigkeit_m"]
C["layers"]["caprock_bottom_thickness_m"] = _d["maechtigkeit_m"]
for _key in ("caprock_top", "caprock_bottom"):
    cr = C["materials"][_key]
    cr["permeability_m2"] = _d["permeabilitaet_m2"]
    cr["porosity"] = _d["porositaet"]
    cr["rho_s_kg_m3"] = _d["dichte_korn_kg_m3"]
    cr["cp_s_J_kgK"] = _d["waermekapazitaet_korn_J_kgK"]
    cr["lambda_s_W_mK"] = _d["waermeleitfaehigkeit_korn_W_mK"]

C["domain"]["size_x_m"] = 2000.0
C["domain"]["size_y_m"] = 400.0
C["wells"]["hot_well_xy"] = (-700.0, 0.0)

_s = FALL["filter_kantenlaenge_m"]
C["wells"]["screen_dx_m"] = _s
C["wells"]["screen_dy_m"] = _s
C["wells"]["screen_top_offset_m"] = 0.0
C["wells"]["screen_bottom_offset_m"] = 0.0
C["wells"]["screen_permeability_m2"] = _k(_a["kf_m_s"])
C["wells"]["production_control"] = "fixed"

P = FALL["monatsleistung_W"]
C["cycles"]["monthly_power_W"] = P
C["cycles"]["monthly_T_inj_K"] = None
C["cycles"]["n_cycles"] = FALL["betriebsjahre"]
C["cycles"]["ramp_days"] = 3.0
DT_REF = FALL["T_injektion_C"] - FALL["T_aquifer_C"]
C["operation"]["mass_flow_rate_kg_s"] = (
    max(abs(p) for p in P) / (C["fluid"]["cp_J_kgK"] * DT_REF))

C["time"]["gravity"] = True
C["dispersion"]["alpha_L_m"] = 2.0
C["dispersion"]["alpha_T_m"] = 0.2

C["fluid"]["beta_1_per_K"] = -4.0e-4

C["operation"]["solid_storage_1_per_Pa"] = 1.0e-9

C["mesh"] = {
    "size_in_well_m":          1.0,
    "size_near_wells_m":       4.0,
    "size_far_m":             35.0,
    "well_size_radius_m":     35.0,
    "well_size_radius_far_m": 160.0,
}

C["solver"]["nonlinear_iter"] = 50

C["time"]["dt_seconds"] = 86400.0
C["time"]["output_every_n_steps"] = 5
C["output"]["out_dir"] = str(Path(__file__).resolve().parent / "ergebnisse_3d")
C["output"]["prefix"] = "ates3d_gw"

if __name__ == "__main__":
    kf = _a["kf_m_s"]
    n = _a["porositaet"]
    i = FALL["gw_gradient"]
    cf = C["fluid"]["cp_J_kgK"]
    R = ((n * RHO_F * cf + (1 - n) * _a["dichte_korn_kg_m3"]
          * _a["waermekapazitaet_korn_J_kgK"]) / (n * RHO_F * cf))
    vt = kf * i / n / R
    print("FALL 2 - ATES 3D MIT regionaler Grundwasserstroemung")
    print(f"  i = {i}, kf = {kf:.2e} m/s  (k = {_k(kf):.3e} m2), R = {R:.2f}")
    print(f"  v_Darcy {kf*i*86400:.2f} m/d | v_Poren {kf*i/n*86400:.2f} m/d "
          f"| v_thermisch {vt*86400:.2f} m/d = {vt*86400*365:.0f} m/a")
    drift = vt * 86400 * 365 * C["cycles"]["n_cycles"]
    platz = C["domain"]["size_x_m"] / 2 - C["wells"]["hot_well_xy"][0]
    print(f"  In {C['cycles']['n_cycles']} Jahren driftet die Front rund "
          f"{drift:.0f} m - Platz stromab: {platz:.0f} m "
          f"({100*drift/platz:.0f} % ausgenutzt)")
    print()
    sys.exit(main())

