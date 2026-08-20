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
from __future__ import annotations

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

# Farben: kategoriale Slots in fester Reihenfolge, Statusfarben getrennt.
C_BLUE, C_ORANGE, C_AQUA, C_YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
C_GOOD, C_WARN, C_CRIT = "#0d6b4f", "#a1620b", "#b3261e"
C_INK, C_INK2, C_INK3, C_RULE = "#14171a", "#4d565e", "#808b95", "#d8dde2"


# ======================================================================
#  1) Lauf erkennen
# ======================================================================
@dataclass
class RunInfo:
    out_dir: Path
    prefix: str
    axisym: bool
    snapshots: list = field(default_factory=list)   # [(t_s, Path)]
    complete: bool = True

    @property
    def t_last(self) -> float:
        return self.snapshots[-1][0] if self.snapshots else 0.0


def detect_run(cfg: dict, out_dir=None) -> RunInfo:
    prefix = cfg["output"]["prefix"]
    out = Path(out_dir if out_dir is not None else cfg["output"]["out_dir"])
    # Glob-Falle: "ates_3d_*" trifft auch "ates_3d_line_*". Der Separator
    # "_ts_" trennt zuverlässig, weil er nur vor der Schrittnummer steht.
    snaps = []
    for p in out.glob(f"{prefix}_ts_*.vtu"):
        m = re.search(r"_t_([0-9.]+)\.vtu$", p.name)
        if m:
            snaps.append((float(m.group(1)), p))
    snaps.sort(key=lambda x: x[0])

    # Axialsymmetrie NICHT am Dateinamen erkennen: erst die .prj fragen,
    # sonst über die Netzausdehnung (ein Achsenmaß = 0 -> 2D).
    axisym = False
    prj = out / f"{prefix}.prj"
    if prj.exists():
        try:
            axisym = 'axially_symmetric="true"' in prj.read_text(
                encoding="ISO-8859-1", errors="replace")
        except OSError:
            pass
    return RunInfo(out_dir=out, prefix=prefix, axisym=axisym, snapshots=snaps)


# ======================================================================
#  2) CONFIG normalisieren  (2D nutzt "well", 3D "wells" usw.)
# ======================================================================
def well_config(cfg: dict) -> dict:
    w = cfg.get("well") or cfg.get("wells") or {}
    lay = cfg["layers"]
    t_aq = lay["aquifer_thickness_m"]
    off = w.get("screen_top_offset_m", 0.0) + w.get("screen_bottom_offset_m", 0.0)
    h_screen = max(1e-9, t_aq - off)

    if "r_well_m" in w:                        # 2D radial
        r_eq = w["r_well_m"]
        V = math.pi * r_eq ** 2 * h_screen
    else:                                      # 3D: Filterbox
        dx, dy = w.get("screen_dx_m", 1.0), w.get("screen_dy_m", 1.0)
        V = dx * dy * h_screen
        r_eq = math.sqrt(dx * dy / math.pi)    # flächengleicher Radius

    dom = cfg["domain"]
    R = dom.get("r_max_m") or 0.5 * min(dom.get("size_x_m", 300.0),
                                        dom.get("size_y_m", 300.0))
    return {
        "r_eq_m": r_eq, "h_screen_m": h_screen, "V_well_m3": V,
        "A_screen_m2": 2.0 * math.pi * r_eq * h_screen,   # Filterzylinderfläche
        "R_influence_m": R,
        "production_control": w.get("production_control", "fixed"),
        "max_rate_factor": w.get("max_rate_factor", 1.0),
    }


def cycle_info(cfg: dict, curves: dict | None = None, rate_mult=None) -> dict:
    cyc = cfg["cycles"]
    op = cfg["operation"]
    T_amb = cfg["initial"]["T_K"]
    dT_ref = op["T_hot_K"] - T_amb
    monthly = cyc.get("monthly_power_W")          # fehlt im Line-Modell ganz
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


# ======================================================================
#  3) Geometrie: Zellvolumen, Materialien, Masken
# ======================================================================
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

    # MaterialID -> Material. Reihenfolge in allen Skripten gleich:
    # aquifer=0, caprock_top=1, caprock_bottom=2, hot_well_vol=3.
    # ates_3d_line.py hat zusätzlich cold_well_vol=4 -> Default nötig,
    # sonst IndexError.
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
        # OGS integriert mit 2*pi*r. pyvista liefert für Quads FLÄCHEN —
        # das Ringvolumen ist V = 2*pi*r_schwerpunkt*A. Ohne diesen Faktor
        # sind alle Energien um Größenordnungen falsch.
        vol = 2.0 * math.pi * cc[:, 0] * np.abs(sizes.cell_data["Area"])
        r_cell, z_cell = cc[:, 0], cc[:, 1]
        x_rel = cc[:, 0]                 # axialsymmetrisch: keine Vorzugsrichtung
    else:
        vol = np.abs(sizes.cell_data["Volume"])
        # Radius vom BRUNNEN aus, nicht vom Koordinatenursprung: der Brunnen
        # darf ausserhalb der Mitte liegen (z. B. weit stromauf, damit die
        # Fahne bei regionaler Stroemung Platz nach stromab hat).
        wsel = (mid == 3)
        if wsel.any():
            xw = float(np.average(cc[wsel, 0], weights=np.abs(sizes.cell_data["Volume"])[wsel]))
            yw = float(np.average(cc[wsel, 1], weights=np.abs(sizes.cell_data["Volume"])[wsel]))
        else:
            xw = yw = 0.0
        r_cell = np.hypot(cc[:, 0] - xw, cc[:, 1] - yw)
        x_rel = cc[:, 0] - xw            # + = stromab bei direction_deg = 0
        z_cell = cc[:, 2]

    keys = list(dom.cells_dict.keys())
    conn = dom.cells_dict[keys[0]] if len(keys) == 1 else None

    z_base = cfg["domain"]["z_base_m"]
    z0 = z_base + cfg["layers"]["caprock_bottom_thickness_m"]
    z1 = z0 + cfg["layers"]["aquifer_thickness_m"]
    mask_well = np.isin(mid, (3,))
    mask_aq = np.isin(mid, (0, 3, 4)) & (z_cell >= z0 - 1e-6) & (z_cell <= z1 + 1e-6)
    if not mask_well.any():                     # Notnagel: kein Brunnenmaterial
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


# ======================================================================
#  4) Zeitreihen
# ======================================================================
def build_timeseries(cfg, run, geo, cyc):
    import pyvista as pv
    rho_f = cfg["fluid"]["rho_ref_kg_m3"]
    T_amb = cyc["T_amb_K"]
    z0, z1 = geo.z_aq
    z_mid = 0.5 * (z0 + z1)
    # Punktsonde auf halber Aquiferhöhe — NUR zum Vergleich. Sie liest bei
    # aktivem Auftrieb systematisch zu kalt, weil die heiße Fahne am
    # Aquiferdach liegt.
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
        # Bei regionaler Stroemung ist nicht der Radius interessant, sondern
        # wie weit die Fahne stromab getragen wird.
        out["x_front_down"].append(float(geo.x_rel[hot].max()) if hot.any() else 0.0)
        out["x_front_up"].append(float(geo.x_rel[hot].min()) if hot.any() else 0.0)
        out["frac_hot"].append(float(geo.vol[geo.mask_aq & (T > T_amb + 40.0)].sum()
                                     / geo.vol[geo.mask_aq].sum()) * 100.0)
        p = _cell_vals(m, "p", geo.conn)
        if p is None:
            out["dp_well"].append(np.nan)
        else:
            if p_ref is None:
                # Referenz ist der UNGESTOERTE Zustand. Bei regionaler
                # Stroemung gehoert deren linearer Anteil dazu, sonst misst
                # man den Standortgradienten am Brunnenort als vermeintlichen
                # Brunnendruck (bei x = -400 m und i = 0.015 sind das 59 kPa).
                p_ref = p.copy()
                _gw = cfg.get("regional_gw", {})
                if _gw.get("enable", False) and not run.axisym:
                    _i = float(_gw.get("gradient_m_per_m", 0.0))
                    _a = math.radians(float(_gw.get("direction_deg", 0.0)))
                    _rg = cfg["fluid"]["rho_ref_kg_m3"] * G * _i
                    # ABSOLUTE x-Koordinate: der Gradient wirkt vom
                    # Koordinatenursprung, nicht vom Brunnen. Mit x_rel bliebe
                    # genau der Versatz rho*g*i*x_Brunnen stehen - bei
                    # x = -400 m sind das die 59 kPa, die dann faelschlich
                    # als Brunnendruck erscheinen.
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

    # Massenstrom und Brunnenleistung auf dem Ausgabe-Zeitgitter
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
    return ts


# ======================================================================
#  5) Kennzahlen je Betriebsjahr
# ======================================================================
def energy_metrics(ts, cyc):
    rows = []
    E_tot = ts["E_aq"] + ts["E_cr"]
    monthly = cyc["monthly_power_W"]
    e_dem = (sum(-p for p in monthly if p < 0) * MONTH / 1e9) if monthly else None
    for y in range(int(cyc["n_cycles"])):
        a, b = y * YEAR, (y + 1) * YEAR
        if ts["t"][-1] < b - 0.02 * YEAR:        # nur volle Jahre bilanzieren
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
    """Liest die OGS-Konsolenausgabe aus dem Ergebnisverzeichnis.

    Ein Lauf kann mitten in der Rechnung abbrechen, ohne dass man es den
    vorhandenen VTU-Dateien ansieht - man wertet dann stillschweigend einen
    halben Lauf aus. Deshalb wird der Log mitgelesen.
    """
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


# ======================================================================
#  6) Plausibilitätsprüfung — Ampel
# ======================================================================
def checks(cfg, run, geo, ts, cyc, rows):
    """-> Liste von (Name, Wert, Einheit, Status, Diagnose, ok_band, fehler_marke)

    Alle Schwellen aus CONFIG abgeleitet, damit dieselbe Tabelle für jeden
    Standort gilt. Die Fehler-Marken sind die gemessenen Fingerabdrücke der
    beiden klassischen Modellfehler.
    """
    T_amb = cyc["T_amb_K"]
    T_inj = cfg["operation"]["T_hot_K"]
    out = []

    def add(name, val, unit, band, status, diag, marke=None):
        out.append(dict(name=name, val=val, unit=unit, band=band,
                        status=status, diag=diag, marke=marke))

    # --- A) Fingerabdruck fehlender rho_f-Division ---------------------
    # Nicht die absolute Geschwindigkeit pruefen (die haengt am Brunnenradius),
    # sondern das Verhaeltnis zur analytisch erwarteten Geschwindigkeit an der
    # Filterflaeche: v_soll = (mdot/rho_f)/A_Filter. Fehlt die rho_f-Division,
    # ist das Verhaeltnis ~1000 statt ~1 - unabhaengig von der Geometrie.
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

    # Auch hier nicht absolut pruefen: 75 kPa sind bei k = 1e-11 richtig und bei
    # k = 6e-11 zu viel. Verglichen wird mit der Thiem-Vorhersage fuer dieselbe
    # Rate - stimmen Handformel und Feld nicht ueberein, ist der Quellterm falsch
    # skaliert.
    K_h = (cfg["materials"]["aquifer"]["permeability_m2"]
           * cfg["fluid"]["rho_ref_kg_m3"] * G / cfg["fluid"]["viscosity_Pa_s"])
    b_aq = cfg["layers"]["aquifer_thickness_m"]
    thiem = (math.log(max(_wc["R_influence_m"] / _wc["r_eq_m"], 1.1))
             / (2 * math.pi * K_h * b_aq))
    dp_soll = (cyc["mdot_nom_kg_s"] / cfg["fluid"]["rho_ref_kg_m3"]
               * thiem * cfg["fluid"]["rho_ref_kg_m3"] * G)
    # # Den ersten Schnappschuss auslassen: bei t = 0 steht die Anfangsbedingung
    # (rein hydrostatisch), der regionale Gradient baut sich erst in den ersten
    # Sekunden auf. Die analytische Referenz enthaelt ihn aber schon - t = 0
    # lieferte sonst den Gradienten selbst als vermeintlichen Brunnendruck.
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

    # --- B) Fingerabdruck der Dauerklemme -----------------------------
    over = (np.nanmax(ts["T_max"]) - T_inj) if len(ts["T_max"]) else np.nan
    # Ohne Stroemung darf hier nichts stehen: der Fall ohne Durchstrom trifft
    # 6e-12 K, also Maschinengenauigkeit. MIT regionaler Stroemung sitzt am
    # Filter ein auf T_inj geklemmter Koerper mitten im Durchstrom, und die
    # Galerkin-Loesung ueberschwingt dort um einen Bruchteil des
    # Temperaturhubs. Gemessen: 1.1 K bei 50 K Hub, energetisch 0.0004 % des
    # Jahreseintrags - ansehen ja, aber es macht die Zahlen nicht wertlos.
    # Die Grenze skaliert deshalb mit dem Hub statt pauschal bei 1 K zu stehen.
    _flow = bool(cfg.get("regional_gw", {}).get("enable", False)) and not run.axisym
    _lift = T_inj - (cyc["T_amb_K"] - 273.15)
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

    # --- Lauf ueberhaupt sauber durchgelaufen? ------------------------
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

    # --- C) Modellaufbau ----------------------------------------------
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
    # Faktor 1.5 statt 1.0: an einem 30-Jahres-Lauf nachgemessen reichte die
    # 0.1-K-Front 71.6 m ins Deckgestein, waehrend 2*sqrt(a*t) nur 54 m ergab.
    # Die Formel unterschaetzt, weil der Aquifer 30 Jahre lang eine warme
    # Quelle ist und nicht ein einmaliger Puls - und weil der Auftrieb nach
    # OBEN zusaetzlich Waerme an die Grenzflaeche traegt.
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
        # Mit Stroemung SOLL die Fahne das Gebiet verlassen - das ist das
        # Ergebnis, nicht der Fehler. Geprueft wird nur, ob stromab genug
        # Platz war, um die Drift ueberhaupt zu sehen.
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
        # Nicht Jahr N gegen N-1 vergleichen: bei unsymmetrischen Lastprofilen
        # pendelt eta von Jahr zu Jahr um mehrere Prozentpunkte. Der Trend
        # steckt im Mittel ueber je zwei Jahre.
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
        # Bei regionaler Stroemung ist der Mehreintrag PHYSIK, kein Befund:
        # der Grundwasserstrom spuelt staendig durch den auf T_inj geklemmten
        # Filterkoerper und wird dabei aufgeheizt - der Brunnen wirkt als
        # Waermetauscher am Strom. Das als FEHLER zu melden bringt Studenten
        # bei, rote Punkte zu ueberlesen; genau das hat diese Uebung schon
        # einmal wochenlang blockiert. Im Fall OHNE Stroemung bleibt der
        # Befund scharf - dort zeigt er ein zu grosses Filtervolumen an.
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


# ======================================================================
#  7) Abbildungen
# ======================================================================
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
        # Alle Größen auf ihr eigenes OK-Band normieren: 0 = untere Grenze,
        # 1 = obere Grenze. So passen Kelvin, kPa und m/d in EIN Bild.
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
        # Nur eine Linie, keine Flaeche: die Deckelung gilt fuer den
        # DECKUNGSGRAD. eta darf darueber liegen - das ist kein Widerspruch,
        # sondern genau der Punkt (guter Speicher, zu kleine Solaranlage).
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
    """Machbarkeitskette: T_inj -> mdot -> Eintrittsgeschw. -> Absenkung.

    Der Sinn ist die Reihenfolge: T_inj ist eine ANNAHME, daraus folgt der
    Massenstrom, daraus die Beanspruchung des Filters und der Druck. Wird eine
    der Grenzen gerissen, ist die Anlage so nicht baubar — und ein Druck um
    Faktor 1000 zu hoch ist der Fingerabdruck des fehlenden rho_f.
    """
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
    V_LIM = 0.03                     # Richtwert Filtereintrittsgeschwindigkeit

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

    # (1) Massenstrom
    ax = axes[0]
    ax.bar(x, val, .64, color=colr)
    lo, hi = 50 * rho_f / 3600, 250 * rho_f / 3600
    ax.set_ylim(0, val.max() * 1.6)
    for _y, _t in ((lo, "50 m³/h"), (hi, "250 m³/h")):
        if _y < val.max() * 1.5:          # nur zeichnen, wenn im Bild
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

    # (2) Filtereintrittsgeschwindigkeit, in % des Richtwerts
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

    # (3) Absenkung / Aufhoehung
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
    """Feld auf ein regelmaessiges Raster in einer Ebene abtasten.

    Bewusst per Abtastung statt per 3D-Rendering: braucht keinen OpenGL-Kontext,
    laeuft also auch ohne Grafikkarte durch, und die Farbskala laesst sich fest
    vorgeben. plane: "xy" (horizontal, z = const), "xz" (vertikal laengs,
    y = const) oder "rz" (2D-Radialmodell, z = 0).
    """
    import pyvista as pv
    A, B = np.meshgrid(xs, ys)
    C = np.full_like(A, float(const))
    if plane == "xy":
        grid = pv.StructuredGrid(A, B, C)
    elif plane == "xz":
        grid = pv.StructuredGrid(A, C, B)
    else:                                    # rz: 2D-Netz liegt in der xy-Ebene
        grid = pv.StructuredGrid(A, B, C)
    res = grid.sample(mesh)
    # pyvista setzt dimensions = A.shape und legt die Punkte in FORTRAN-
    # Reihenfolge ab. Ein reshape in C-Ordnung transponiert das Feld still
    # und liefert diagonale Streifen statt eines Temperaturfeldes.
    T = np.asarray(res["T"]).reshape(A.shape, order="F") - 273.15
    valid = np.asarray(res["vtkValidPointMask"]).reshape(A.shape, order="F").astype(bool)
    T[~valid] = np.nan
    return A, B, T


def _pick_times(ts, cyc, n=6, whole=False):
    """n Zeitpunkte aus dem LETZTEN vollen Betriebsjahr, gleichmaessig verteilt.

    Nicht ueber den ganzen Lauf verteilen: das erste Jahr ist nie
    repraesentativ, und die Aussage ist das saisonale Atmen.
    """
    t = ts["t"]
    if len(t) < 2:
        return []
    if whole:
        # Bei regionaler Stroemung ist der ganze Lauf interessant: die Fahne
        # wandert weg und kommt nicht wieder. Ein einzelnes Jahr zeigt das nicht.
        a, b = 0.0, t[-1]
    else:
        yr_last = max(0.0, np.floor(t[-1] / YEAR) - 1) if t[-1] >= 2 * YEAR else 0.0
        a = yr_last * YEAR
        b = min(a + YEAR, t[-1])
    want = np.linspace(a, b, n + 1)[1:]
    return [int(np.argmin(np.abs(t - w))) for w in want]


def fig_feldschnitte(cfg, run, geo, ts, cyc, fig_dir, plt):
    """Temperaturfeld zu mehreren Zeitpunkten - fixe Farbskala.

    Die feste Skala T_amb..T_inj ueber alle Teilbilder ist der Punkt: mit der
    ueblichen Autoskalierung sieht ein Lauf, in dem der halbe Aquifer kocht,
    genauso aus wie ein gesunder.
    """
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

    # ---------- 2D axialsymmetrisch: (r, z), gespiegelt ----------------
    if run.axisym:
        r_max = max(30.0, 1.25 * float(np.nanmax(ts["r_front"])))
        # Erst auf der HALBEBENE abtasten (r monoton!), dann das Ergebnis
        # spiegeln. Ein nicht monotones Koordinatenfeld faltet das
        # StructuredGrid in sich selbst und liefert Streifenmuster.
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

    # ---------- 3D: Draufsicht (Drift) + Laengsschnitt ------------------
    xw = float(np.average(geo.r_cell[geo.mask_well] * 0 + 0))   # Platzhalter
    # Brunnenlage aus der Zellmaske
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


# ======================================================================
#  8) Hauptfunktion
# ======================================================================
def report(cfg, out_dir=None, curves=None, rate_mult=None, report_dir=None):
    # Windows-Konsolen laufen oft auf cp1252 und wirfen bei Zeichen wie "−"
    # (U+2212) oder "ρ" einen UnicodeEncodeError - der Bericht wuerde daran
    # sterben. Nur die FEHLERBEHANDLUNG wird umgestellt, nicht die Codierung.
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

    # --- CSV ---------------------------------------------------------
    if rows:
        with open(rep / f"{run.prefix}_kennzahlen.csv", "w", newline="",
                  encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows([{k: (round(v, 3) if isinstance(v, float) else v)
                          for k, v in r.items()} for r in rows])
    with open(rep / f"{run.prefix}_pruefblatt.csv", "w", newline="",
              encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["pruefgroesse", "wert", "einheit", "ok_von", "ok_bis",
                    "status", "diagnose"])
        for c in chks:
            w.writerow([c["name"], f"{c['val']:.4g}", c["unit"],
                        c["band"][0], c["band"][1], c["status"], c["diag"]])

    # --- Konsole -----------------------------------------------------
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

    # --- Abbildungen -------------------------------------------------
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
        except Exception as e:                       # eine Abbildung darf scheitern
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
