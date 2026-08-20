#!/usr/bin/env python3
"""ATES 2D radialsymmetrisch, OHNE regionale Grundwasserstroemung.

Die Frage: Was leistet der Aquiferspeicher am Standort im guenstigsten Fall?
Ohne Hintergrundstroemung bleibt die eingespeicherte Waerme am Brunnen; das ist
die Obergrenze dessen, was der Speicher liefern kann.

Ein axialsymmetrisches Modell KANN keine gerichtete Stroemung abbilden - das
ist hier kein Mangel, sondern genau der gewuenschte Fall. Die Gegenprobe mit
Stroemung ist ates_3D.py.

    python ates_2D.py               # 30 Betriebsjahre (rund 1 h 45 min)
    python ates_2D.py --years 5     # kurzer Durchlauf zum Ausprobieren
    python ates_2D.py --no-run      # nur Netz + .prj ansehen, nicht rechnen

Ergebnisse landen neben dieser Datei in ergebnisse_2d/:
    ergebnisse_2d/ates2d.pvd       fuer ParaView (zieht alle Zeitschritte)
    ergebnisse_2d/figures/         Abbildungen und Pruefblatt
    ergebnisse_2d/*_kennzahlen.csv Rueckgewinnungsgrad, Deckungsgrad je Jahr

ZUM BEARBEITEN gibt es genau EINEN Block: das Dict FALL weiter unten.
Netz, Zeitschritt, Loeser und Ausgabe stehen darunter, sind eingestellt und
begruendet - die muss man fuer diese Uebung nicht anfassen. Der Motor liegt
im Unterordner modell/ und wird ebenfalls nicht angefasst.
"""
from __future__ import annotations

# ######################################################################
#
#   H I E R   S C H R A U B E N
#
#   Standort und Lastfall. Alles andere ist fertig eingestellt.
#
# ######################################################################

FALL = {

    # ------------------------------------------------------------------
    #  Lastprofil: zwoelf Monatsleistungen in Watt
    # ------------------------------------------------------------------
    #  P > 0 = einspeichern (laden),  P < 0 = foerdern (entladen).
    #
    #  Das ist die wichtigste Eingabe der ganzen Rechnung. Die Summe der
    #  Einspeisung geteilt durch die Summe der Entnahme DECKELT den
    #  Deckungsgrad, bevor eine einzige Zelle gerechnet wird:
    #  bei diesem Profil 3061/6721 = 45.5 %. Mehr als 45.5 % des Bedarfs
    #  kann der Speicher auch bei einem Rueckgewinnungsgrad von 100 %
    #  nicht liefern. Das Pruefblatt weist die Deckelung aus.
    #
    #  Die Werte unten sind Solarueberschuss minus Waermebedarf aus der
    #  Solarthermie-Uebung (halbe Krankenhauslast).
    "monatsleistung_W": [
        -595_850.0,   # Jan
        -523_490.0,   # Feb
        -336_450.0,   # Mrz
          +4_975.0,   # Apr
        +162_190.0,   # Mai
        +343_595.0,   # Jun
        +346_995.0,   # Jul
        +282_210.0,   # Aug
         +24_135.0,   # Sep
        -281_920.0,   # Okt
        -370_970.0,   # Nov
        -447_005.0,   # Dez
    ],

    # ------------------------------------------------------------------
    #  Betrieb
    # ------------------------------------------------------------------
    "T_injektion_C":   60.0,   # womit geladen wird
    "T_aquifer_C":     10.0,   # ungestoerte Untergrundtemperatur
    "betriebsjahre":     30,   # ueber --years ueberschreibbar
    #
    #  Einen Knopf fuer die Auslegungsspreizung gibt es bewusst NICHT.
    #  Das Modell rechnet den Massenstrom im Monatsprofil-Modus immer als
    #      mdot = P_max / (c_f * (T_injektion - T_aquifer))
    #  und liest operation.mass_flow_rate_kg_s dabei gar nicht. Dahinter
    #  steckt die Annahme, dass das gefoerderte Wasser bei
    #  Aquifertemperatur zurueckgespeist wird. Eine andere Spreizung ist
    #  mit diesem Skript nicht rechenbar - ein Knopf dafuer waere tot und
    #  wuerde nur die Konsolenausgabe veraendern.

    # ------------------------------------------------------------------
    #  Aquifer  (Standort aus der Abgabe)
    # ------------------------------------------------------------------
    "aquifer": {
        "maechtigkeit_m":                 38.0,
        # Durchlaessigkeitsbeiwert, der gemessene Standortwert. Die
        # Umrechnung auf die Permeabilitaet k = kf*mu/(rho*g) macht das
        # Skript weiter unten selbst - und zwar mit der Viskositaet, die
        # das Modell bei Aquifertemperatur wirklich verwendet. Wer hier
        # k statt kf eintraegt, baut sich denselben Fehler ein, der in
        # der Vorgaengerfassung steckte (siehe Kommentar bei _k()).
        "kf_m_s":                       6.0e-4,
        "porositaet":                   0.1191,
        "dichte_korn_kg_m3":            2760.0,
        "waermekapazitaet_korn_J_kgK":   793.0,
        "waermeleitfaehigkeit_korn_W_mK": 2.28,
    },

    # ------------------------------------------------------------------
    #  Deckgestein  (oben und unten gleich aufgebaut)
    # ------------------------------------------------------------------
    "deckgestein": {
        # 120 m je Seite. Gemessen an einem Vorlauf mit 80 m: die
        # Waermeleitfront (dT > 0.1 K) stand nach 30 Jahren 71.6 m weit
        # im oberen Deckgestein - nur 8.4 m vor dem auf T_amb fixierten
        # Modellrand. Nach unten 52.9 m. Die Asymmetrie ist der Auftrieb.
        "maechtigkeit_m":               120.0,
        "permeabilitaet_m2":          2.1e-16,
        # k_vertikal / k_horizontal. ACHTUNG beim Vergleich mit Fall 2:
        # dieser Knopf existiert nur im 2D-Modell. Das 3D-Modell schreibt
        # eine skalare Permeabilitaet, sein Deckgestein ist zwingend isotrop
        # und vertikal damit zehnmal durchlaessiger als hier. Wer diese Zeile
        # ins 3D-Handout kopiert, schreibt sie in ein Dict, das nie gelesen
        # wird - kein Fehler, keine Warnung, keine Wirkung.
        "anisotropie_vertikal":           0.1,
        "porositaet":                    0.05,
        "dichte_korn_kg_m3":           2700.0,
        "waermekapazitaet_korn_J_kgK":  900.0,
        "waermeleitfaehigkeit_korn_W_mK": 2.0,
    },

    # ------------------------------------------------------------------
    #  Brunnen
    # ------------------------------------------------------------------
    "brunnenradius_m": 0.1,
    # Die Filterstrecke geht ueber die volle Aquifermaechtigkeit. Anders
    # kann das 2D-Netz es nicht: es hat nur eine Zellspalte am Brunnen.

}

# ####################################################################
#
#   A B   H I E R   S T E H T   D E R   R E C H E N K E R N
#   (Netz, OGS-Projektdatei, Loesen, Pruefbericht) - nicht anfassen.
#
# ####################################################################
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

# ####################################################################
#  Der Pruefbericht steht oben in dieser Datei. Der Rechenkern
#  importiert ihn als Modul "ates_report" - also melden wir eines an,
#  statt den Kern anzufassen.
# ####################################################################
import sys as _sys
import types as _types
_mod = _types.ModuleType("ates_report")
for _n in ("auto_report", "report"):
    setattr(_mod, _n, globals()[_n])
_sys.modules["ates_report"] = _mod
#!/usr/bin/env python3
"""
ATES 2D radialsymmetrisch — Einführungsbeispiel für OpenGeoSys 6.

Ein einzelner Brunnen sitzt auf der Symmetrieachse (r = 0). Das 3D-Problem
wird unter Annahme rotations­symmetrischer Lösung auf eine 2D-Aufgabe in
der (r, z)-Halbebene reduziert. OGS-6 berücksichtigt die Achssymmetrie
automatisch, wenn `axially_symmetric="true"` im Mesh-Eintrag steht.

VORTEIL: Bei gleicher Auflösung 100-1000x weniger Zellen als 3D
         -> Sekunden bis wenige Minuten Laufzeit.

VERWENDUNG
----------
    python ates_radial_2d.py             # Mesh + Sim + Plots
    python ates_radial_2d.py --no-run    # nur Setup
    python ates_radial_2d.py --no-mesh   # nur .prj
    python ates_radial_2d.py --no-plots  # ohne Auto-Plots

Konventionen
------------
    x = r  (Radialkoordinate, x >= 0)
    y = z  (Vertikalkoordinate)
"""

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

import gmsh
import numpy as np

def _safe_name(name):
    """Physical-Group-Namen plattform- und dateisystemsicher machen.

    Namen koennen Zeichen enthalten, die in Dateinamen unzulaessig sind
    (Slash "/", Backslash, Leerzeichen, Umlaute, ...). Ein Slash wuerde
    von pyvista als Ordnertrennung gedeutet -> FileNotFoundError. Nur
    ASCII-Buchstaben/Ziffern sowie . _ - bleiben erhalten; alles andere
    wird durch "_" ersetzt. So laeuft es auf Windows, Linux und macOS
    gleichermassen.
    """
    keep = ("abcdefghijklmnopqrstuvwxyz"
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")
    return "".join(c if c in keep else "_" for c in str(name))


# --- ogstools >=0.8 Kompatibilitäts-Shim für die alte msh2vtu-API ---
def msh2vtu(filename, output_path, output_prefix, dim, reindex=True, log_level="WARNING"):
    import ogstools as ot
    from pathlib import Path as _P
    meshes = ot.Meshes.from_gmsh(
        filename=str(filename), dim=dim, reindex=reindex, log=False
    )
    output_path = _P(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    for name, mesh in meshes.items():
        if name == "domain":
            fname = f"{output_prefix}_domain.vtu"
        else:
            fname = f"{output_prefix}_physical_group_{_safe_name(name)}.vtu"
        mesh.save(str(output_path / fname), binary=True)


# ======================================================================
#  CONFIG  --  alle einstellbaren Parameter
# ======================================================================
CONFIG: dict = {

    # --- Domäne ---------------------------------------------------------
    # An MOOSE-ATES angelehnt (mooseframework.inl.gov/.../ates.html):
    #   Modellradius 1000 m, Aquifer 20 m, Kaprock 2×40 m, Brunnen 0,1 m.
    "domain": {
        "r_max_m":  1000.0,
        "z_base_m":    0.0,
    },

    # --- Schichtdicken --------------------------------------------------
    # Deckgestein dick genug, dass die 30-Jahre-Wärmeleitung (~35 m) den auf T_amb
    # fixierten Außenrand nicht erreicht -> 80 m je Seite (klare Reserve).
    "layers": {
        "caprock_bottom_thickness_m": 80.0,
        "aquifer_thickness_m":        20.0,
        "caprock_top_thickness_m":    80.0,
    },

    # --- Brunnen --------------------------------------------------------
    "well": {
        "r_well_m":                0.1,    # Brunnenradius [m] (MOOSE: bh_r=0.1)
        "screen_top_offset_m":     0.0,    # Filter über volle Aquiferhöhe
        "screen_bottom_offset_m":  0.0,
        "screen_permeability_m2":  1.0e-11,# = Aquifer (Filter voll durchlässig)
        # Förderregelung:
        #   "fixed"  : Pumprate = |P|/(c_p·ΔT_ref) fest (wie MOOSE); gelieferte
        #              Leistung sinkt mit der Fördertemperatur.
        #   "demand" : bedarfsgeführt — Rate wird (monatlich, iterativ) so
        #              nachgeregelt, dass die geforderte Monatsleistung gehalten
        #              wird: ṁ = |P|/(c_p·(T̄_prod−T_amb)), gedeckelt auf
        #              max_rate_factor × Nennrate.
        "production_control":      "demand",
        "max_rate_factor":         6.0,
        "demand_iterations":       5,      # Fixpunkt-Iterationen der Bedarfsführung
    },

    # --- Mesh-Auflösung -------------------------------------------------
    # Strukturiertes Schichtnetz (Quads, transfinit). Zellen sind an den
    # Aquifer-/Deckgestein-Grenzflächen konform -> scharfe Schichtgrenzen.
    #   radial:  Brunnen | fein bis r_fine | grob bis r_max (geometrisch)
    #   vertikal: Aquifer gleichmäßig; Deckgestein fein an der Grenzfläche,
    #             gröber zum Außenrand (geometrischer Bias).
    "mesh": {
        "r_fine_m":            70.0,   # radial fein bis hierher (deckt Fahne + Rand ab)
        "n_r_well":             4,     # Zellen über den Brunnenradius (0.1 m)
        "n_r_fine":            84,     # radiale Zellen [r_well .. r_fine]  (~0.1→2 m)
        "n_r_far":             34,     # radiale Zellen [r_fine .. r_max]   (plumefern)
        "bias_r_fine":         1.045,  # radiale Progression (fein am Brunnen)
        "bias_r_far":          1.05,   # radiale Progression (fein bei r_fine)
        "n_z_aquifer":         28,     # vertikale Zellen im Aquifer (~0.7 m)
        # Deckgestein: feine Bandzone an der Grenzfläche (löst die Wärmeleitfront
        # auf) + grobe Zone zum Außenrand. Tiefe der feinen Zone konfigurierbar.
        "caprock_fine_depth_m": 45.0,  # feine Zone reicht so tief ins Deckgestein
        "n_z_caprock_fine":     30,    # Zellen in der feinen Deckgestein-Zone (~1.5 m)
        "n_z_caprock_coarse":   10,    # Zellen in der groben Zone (biased)
        "bias_caprock":         1.28,  # Progression der groben Zone (fein zur Front)
    },

    # --- Materialien (MOOSE-Werte) -------------------------------------
    # Anisotrope Permeabilität: (k_hor, k_ver). Aquifer & Kaprock gleiche
    # thermische Eigenschaften, nur Permeabilität unterscheidet sich.
    "materials": {
        "aquifer": {
            "permeability_m2":     1.0e-11,   # horizontal
            "permeability_ver_m2": 2.0e-12,   # vertikal (MOOSE: 5× kleiner)
            "porosity":        0.25,
            "rho_s_kg_m3":  2650.0,
            "cp_s_J_kgK":    800.0,
            "lambda_s_W_mK":   3.0,
        },
        "caprock_top": {
            "permeability_m2":     1.0e-16,
            "permeability_ver_m2": 1.0e-17,
            "porosity":        0.25,
            "rho_s_kg_m3":  2650.0,
            "cp_s_J_kgK":    800.0,
            "lambda_s_W_mK":   3.0,
        },
        "caprock_bottom": {
            "permeability_m2":     1.0e-16,
            "permeability_ver_m2": 1.0e-17,
            "porosity":        0.25,
            "rho_s_kg_m3":  2650.0,
            "cp_s_J_kgK":    800.0,
            "lambda_s_W_mK":   3.0,
        },
    },

    # --- Fluid ----------------------------------------------------------
    # Für Auftrieb: T-abhängige Dichte ρ(T)=ρ_ref·(1+β·(T−T_ref)) mit β<0.
    # (Boussinesq-artig; ersetzt MOOSEs Water97-Tabelle.)
    "fluid": {
        "rho_ref_kg_m3":  1000.0,
        "T_ref_K":         283.15,
        "beta_1_per_K":   -4.0e-4,   # thermische Ausdehnung (heißes Wasser leichter)
        # Viskosität T-abhängig μ(T)=μ_ref·(1+slope·(T−T_ref)); verstärkt den
        # Auftrieb (heißes Wasser dünnflüssiger). μ_ref bei T_ref=10 °C ≈ 1,3e-3,
        # slope so, dass μ(60 °C) ≈ 0,47e-3 (wie echtes Wasser). None → konstant.
        "viscosity_Pa_s":       1.3e-3,
        "visc_slope_1_per_K":  -1.28e-2,
        "cp_J_kgK":        4180.0,
        "lambda_W_mK":     0.6,
    },

    # --- Dispersivität --------------------------------------------------
    "dispersion": {
        "alpha_L_m": 5.0,
        "alpha_T_m": 0.5,
    },

    # --- Anfangsbedingungen ---------------------------------------------
    # T_amb = Umgebungstemperatur des Aquifers (= fluid.T_ref_K, damit
    # ρ(T_amb)=ρ_ref). p_top_Pa = Druck an der Modelloberkante; darunter
    # hydrostatisch (bei aktivierter Gravitation, siehe build_prj).
    "initial": {
        "T_K":  283.15,      # 10 °C Umgebung
        "p_Pa": 1.0e5,       # Referenzdruck Oberkante (hydrostatisch nach unten)
        "T_surface_K":                  283.15,
        "geothermal_gradient_K_per_m":  0.0,
    },

    # --- Betrieb --------------------------------------------------------
    # T_hot_K = Injektionstemperatur T_inj. Der Massenstrom folgt aus der
    # Monatsleistung:  ṁ = P / ( c_p,f · (T_inj − T_amb) ).
    "operation": {
        "mass_flow_rate_kg_s": 10.0,       # nur Referenz im 4-Phasen-Modus
        "T_hot_K":  333.15,                # T_inj = 60 °C
        "T_cold_K": 283.15,
        "fluid_storage_1_per_Pa": 4.5e-10,
        "solid_storage_1_per_Pa": 1.0e-10,
    },

    # ------------------------------------------------------------------
    # ZYKLEN
    # ------------------------------------------------------------------
    # Zwei Modi (alternativ):
    #
    # A) 4-Phasen-Zyklus (Default):
    #    Pro Zyklus 4 aufeinander folgende Phasen mit fester Dauer.
    #    Massenstrom = ±operation.mass_flow_rate_kg_s; T_inj = T_hot bzw. T_cold.
    #
    # B) Monatsprofil (überschreibt A, falls aktiviert):
    #    Liste von 12 Monats-Speicher­leistungen [W]. Positiv = laden
    #    (Wärme in den Aquifer), negativ = fördern, 0 = Stillstand.
    #    Daraus wird pro Monat der Massenstrom berechnet:
    #        ṁ_Monat = P_Monat / ( c_p,f · ( T_inj − T_0 ) )
    #    Vorlauf-Temperatur T_inj wahlweise konstant (Default: T_hot_K bei
    #    Beladung, T_cold_K bei Förderung) oder pro Monat über
    #    monthly_T_inj_K = [T_Jan, …, T_Dez].
    #    operation.mass_flow_rate_kg_s wird dann zur Referenz-/Skalierungsgröße
    #    (z. B. Maximalwert eintragen).
    #    n_cycles entspricht der Anzahl Betriebs­jahre.
    # ------------------------------------------------------------------
    "cycles": {
        "n_cycles":                          2,       # Anzahl Zyklen (Modus A) bzw. Betriebsjahre (Modus B)
        # --- Modus A: 4-Phasen-Zyklus ---
        "charge_days":                      91.25,   # Beladung (Tage)
        "storage_after_charge_days":        91.25,   # Pause nach Beladung
        "discharge_days":                   91.25,   # Förderung
        "storage_after_discharge_days":     91.25,   # Pause nach Förderung
        "ramp_days":                         3.0,    # Übergangsrampe zwischen Phasen/Monaten
        # --- Modus B: Monatsprofil (AKTIV; auf None für Modus A) ---
        # Monatliche Speicherleistung P [W] (Jan … Dez).
        #   P > 0  -> Beladung  (Wärme in den Aquifer injizieren)
        #   P < 0  -> Förderung (Wärme entnehmen)
        #   P = 0  -> Stillstand
        # Aus P wird automatisch die Pumprate berechnet:
        #     ṁ_Monat = P_Monat / ( c_p,f · ΔT_ref )   mit ΔT_ref = T_hot_K − T0
        # und der Wärmequellterm  q_h = P_Monat / V_well  [W/m³] eingeprägt.
        # Das hier hinterlegte Saisonprofil ist über das Jahr bilanziert (ΣP = 0).
        # MOOSE-Skala: Spitze 2 MW → ṁ ≈ P/(c_p·ΔT) ≈ 2e6/(4180·50) ≈ 9,6 kg/s
        # (vgl. MOOSE ~12,7 kg/s). ΔT = T_inj − T_amb = 60−10 = 50 K.
        "monthly_power_W": [
            -1_600_000,  # Jan  Förderung (Winter)
            -1_400_000,  # Feb
              -500_000,  # Mär
                     0,  # Apr  Übergang
            +1_000_000,  # Mai  Beladung (Sommer)
            +1_800_000,  # Jun
            +2_000_000,  # Jul
            +1_600_000,  # Aug
              +500_000,  # Sep
              -700_000,  # Okt  Förderung (Herbst/Winter)
            -1_300_000,  # Nov
            -1_400_000,  # Dez
        ],
        # Optional: monatliche Vorlauf­temperatur [K] bei Beladung; None → T_hot_K.
        "monthly_T_inj_K":                  None,
    },

    # --- Zeit + Output --------------------------------------------------
    "time": {
        # Feine, ~tägliche Zeitauflösung: hält die Courant-Zahl am Brunnen O(1)
        # und vermeidet die numerische Dispersion grober (wöchentlicher) Schritte.
        "dt_seconds":            1 * 86400.0,
        # Ausgabe nicht jeden Tag (bei 30 Jahren sonst ~11 000 VTU-Dateien),
        # sondern alle N Rechenschritte -> löst das Monatsprofil weiterhin fein auf.
        "output_every_n_steps":  5,
        "gravity":               True,     # Auftrieb (heißes Wasser steigt) wie MOOSE
    },
    "output": {
        "prefix":    "ates_radial_2d",
        "out_dir":   "out",
        "variables": ["T", "p", "darcy_velocity"],
    },
    "solver": {
        # solver_type: "BiCGSTAB" (iterativ, schnell bei feinem Netz) oder
        # "SparseLU" (direkt, robust). Auf dem strukturierten Quad-Netz ist
        # BiCGSTAB+ILUT gut konditioniert.
        "solver_type":     "BiCGSTAB",
        "precon_type":     "ILUT",
        "linear_tol":      1.0e-10,
        "linear_iter":     10000,
        "nonlinear_iter":  20,
        "rel_tol_T":       1.0e-4,
        "rel_tol_p":       1.0e-4,
    },
}

DAY = 86400.0


# ======================================================================
#  1) Mesh-Generierung (2D radial)
# ======================================================================
def build_mesh(cfg: dict, out_dir: Path) -> Path:
    """Strukturiertes 2D-Schichtnetz (Quads, transfinit) in (r, z)-Halbebene.

    3 vertikale Schichten (Deckgestein unten / Aquifer / Deckgestein oben) ×
    3 radiale Spalten (Brunnen / fein / grob) = 9 transfinite Quad-Flächen.
    Die Zellen sind an den Schichtgrenzen KONFORM -> scharfe Aquifer-/
    Deckgestein-Übergänge. Vertikal fein an den Grenzflächen (geometrischer
    Bias), radial fein am Brunnen. x = r, y = z.
    """
    msh_path = out_dir / f"{cfg['output']['prefix']}.msh"
    r_max  = cfg["domain"]["r_max_m"]
    z_base = cfg["domain"]["z_base_m"]
    t_cb = cfg["layers"]["caprock_bottom_thickness_m"]
    t_aq = cfg["layers"]["aquifer_thickness_m"]
    t_ct = cfg["layers"]["caprock_top_thickness_m"]
    z_aq_bot = z_base + t_cb
    z_aq_top = z_aq_bot + t_aq
    z_top    = z_aq_top + t_ct
    r_well = cfg["well"]["r_well_m"]
    m = cfg["mesh"]

    # --- radiale Segmente: (r0, r1, n_cells, progression) ---
    # Optional vierteilig: Brunnen | Nahfeld (stark gestaffelt) | Fahnenband
    # (nahezu gleichmaessig) | Fernfeld. Der Sinn: mit EINER geometrischen
    # Progression vom Brunnen bis zum Fahnenrand wird entweder der Brunnen zu
    # grob oder der Fahnenrand zu fein. Sobald mesh.r_near_m gesetzt ist, wird
    # das Nahfeld vom Fahnenband getrennt. Ohne den Schluessel bleibt es beim
    # alten Dreisatz (abwaertskompatibel).
    if m.get("r_near_m"):
        rsegs = [(0.0,           r_well,        m["n_r_well"], 1.0),
                 (r_well,        m["r_near_m"], m["n_r_near"],
                  m.get("bias_r_near", 1.05)),
                 (m["r_near_m"], m["r_fine_m"], m["n_r_fine"],
                  m.get("bias_r_fine", 1.006)),
                 (m["r_fine_m"], r_max,         m["n_r_far"],  m["bias_r_far"])]
    else:
        rsegs = [(0.0,      r_well,        m["n_r_well"], 1.0),
                 (r_well,   m["r_fine_m"], m["n_r_fine"], m["bias_r_fine"]),
                 (m["r_fine_m"], r_max,    m["n_r_far"],  m["bias_r_far"])]

    # --- vertikale Segmente: Deckgestein grob | fein | Aquifer | fein | grob ---
    bc = m["bias_caprock"]
    fd_cb = min(m["caprock_fine_depth_m"], t_cb)
    fd_ct = min(m["caprock_fine_depth_m"], t_ct)
    zsegs = []
    if t_cb - fd_cb > 1e-6:      # untere grobe Zone (fein zur Grenzfläche oben -> prog<1)
        zsegs.append((z_base, z_aq_bot - fd_cb, m["n_z_caprock_coarse"], 1.0 / bc))
    zsegs.append((z_aq_bot - fd_cb, z_aq_bot,   m["n_z_caprock_fine"], 1.0))   # untere feine Zone
    aq_j = len(zsegs)
    zsegs.append((z_aq_bot, z_aq_top, m["n_z_aquifer"], 1.0))                  # Aquifer
    zsegs.append((z_aq_top, z_aq_top + fd_ct, m["n_z_caprock_fine"], 1.0))     # obere feine Zone
    if t_ct - fd_ct > 1e-6:      # obere grobe Zone (fein zur Grenzfläche unten -> prog>1)
        zsegs.append((z_aq_top + fd_ct, z_top, m["n_z_caprock_coarse"], bc))

    rx = [rsegs[0][0]] + [s[1] for s in rsegs];  NI = len(rsegs)
    zy = [zsegs[0][0]] + [s[1] for s in zsegs];  NJ = len(zsegs)
    nx = [s[2] + 1 for s in rsegs];  px = [s[3] for s in rsegs]
    ny = [s[2] + 1 for s in zsegs];  py = [s[3] for s in zsegs]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("ates_radial_2d")
    geo = gmsh.model.geo

    # Gitterpunkte / Kanten / Flächen (strukturiertes NI×NJ-Gitter)
    P = {(i, j): geo.addPoint(rx[i], zy[j], 0.0)
         for i in range(NI + 1) for j in range(NJ + 1)}
    hL = {(i, j): geo.addLine(P[(i, j)], P[(i + 1, j)])
          for i in range(NI) for j in range(NJ + 1)}
    vL = {(i, j): geo.addLine(P[(i, j)], P[(i, j + 1)])
          for i in range(NI + 1) for j in range(NJ)}
    S = {}
    for i in range(NI):
        for j in range(NJ):
            cl = geo.addCurveLoop([hL[(i, j)], vL[(i + 1, j)], -hL[(i, j + 1)], -vL[(i, j)]])
            S[(i, j)] = geo.addPlaneSurface([cl])
    geo.synchronize()

    for i in range(NI):
        for j in range(NJ + 1):
            geo.mesh.setTransfiniteCurve(hL[(i, j)], nx[i], "Progression", px[i])
    for i in range(NI + 1):
        for j in range(NJ):
            geo.mesh.setTransfiniteCurve(vL[(i, j)], ny[j], "Progression", py[j])
    for i in range(NI):
        for j in range(NJ):
            geo.mesh.setTransfiniteSurface(S[(i, j)])
            geo.mesh.setRecombine(2, S[(i, j)])
    geo.synchronize()

    # Physical Groups (Reihenfolge -> MaterialID nach reindex: aq=0,ct=1,cb=2,well=3)
    aquifer = [S[(i, aq_j)] for i in range(1, NI)]          # ohne Brunnenspalte (i=0)
    well    = [S[(0, aq_j)]]
    cb      = [S[(i, j)] for i in range(NI) for j in range(aq_j)]        # unter Aquifer
    ct      = [S[(i, j)] for i in range(NI) for j in range(aq_j + 1, NJ)]# über Aquifer
    gmsh.model.addPhysicalGroup(2, aquifer, tag=1, name="aquifer")
    gmsh.model.addPhysicalGroup(2, ct,      tag=2, name="caprock_top")
    gmsh.model.addPhysicalGroup(2, cb,      tag=3, name="caprock_bottom")
    gmsh.model.addPhysicalGroup(2, well,    tag=4, name="hot_well_vol")
    # Randkanten
    gmsh.model.addPhysicalGroup(1, [hL[(i, NJ)] for i in range(NI)], tag=10, name="top")
    gmsh.model.addPhysicalGroup(1, [hL[(i, 0)]  for i in range(NI)], tag=11, name="bottom")
    gmsh.model.addPhysicalGroup(1, [vL[(NI, j)] for j in range(NJ)], tag=12, name="far")
    gmsh.model.addPhysicalGroup(1, [vL[(0, j)]  for j in range(NJ)], tag=13, name="axis")
    gmsh.model.addPhysicalGroup(1, [vL[(1, aq_j)]], tag=14, name="hot_well_surf")

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.model.mesh.generate(2)
    gmsh.write(str(msh_path))
    gmsh.finalize()
    return msh_path


def convert_mesh(cfg: dict, msh_path: Path, out_dir: Path) -> dict[str, str]:
    pass  # msh2vtu provided at module level
    prefix = cfg["output"]["prefix"]
    msh2vtu(filename=msh_path, output_path=out_dir, output_prefix=prefix,
            dim=2, reindex=True, log_level="WARNING")
    return {
        "domain":         f"{prefix}_domain.vtu",
        "top":            f"{prefix}_physical_group_top.vtu",
        "bottom":         f"{prefix}_physical_group_bottom.vtu",
        "far":            f"{prefix}_physical_group_far.vtu",
        "hot_well_vol":   f"{prefix}_physical_group_hot_well_vol.vtu",
        "hot_well_surf":  f"{prefix}_physical_group_hot_well_surf.vtu",
    }


# ======================================================================
#  2) Saisonale Kurven
# ======================================================================
def build_cycle_curves(cfg: dict, rate_mult=None) -> dict:
    """Erzeuge die normierte Pump-/Leistungskurve g(t) ∈ [−1, +1] und die
    Beladungs-Intervalle für die (zeitlich begrenzte) Injektions-Dirichlet-BC.

    rate_mult: optionale Liste von 12 Faktoren; skaliert die Förder-Pumprate je
    Monat (nur P<0). Für die bedarfsgeführte Förderung (well.production_control
    == 'demand') liefert die äußere Iteration diese Faktoren; None → alle 1.

    Brunnenmodell (physikalisch korrekt für OGS-HT, nicht-konservative Form):

      * Massenquelle (p-Gleichung), IMMER aktiv, kurvengesteuert:
            q_m(t) = g(t) · (ṁ_nom / V_well)   [kg/m³/s]
        Kopplung Leistung → Pumprate über feste Referenzspreizung
        ΔT_ref = T_hot − T0:   ṁ_nom = P_nom / (c_p,f · ΔT_ref).
        g > 0 = Injektion (Massenquelle), g < 0 = Förderung (Massensenke).

      * KEINE explizite Wärmequelle. Die Enthalpie wird korrekt über den in
        der nicht-konservativen Advektion enthaltenen Term c_f·T·Q_m
        transportiert:
          - Beladung: Injektionstemperatur wird per Dirichlet T=T_inj am
            Brunnen vorgegeben (nur in den Beladungs-Intervallen aktiv) →
            die Massenquelle trägt Enthalpie bei genau T_inj ein.
          - Förderung: KEIN Dirichlet. Die Massensenke entzieht Enthalpie bei
            der lokalen, dynamisch berechneten Speichertemperatur → die
            Entnahmetemperatur ergibt sich physikalisch von selbst.

    Rückgabe: cycle_power (Massenkurve) und charge_intervals als Liste von
    (t_start, t_end, T_inj_K).
    """
    cyc   = cfg["cycles"]
    n     = cyc["n_cycles"]
    ramp  = max(60.0, cyc["ramp_days"] * DAY)

    T0       = cfg["initial"]["T_K"]
    T_hot    = cfg["operation"]["T_hot_K"]
    cp_f     = cfg["fluid"]["cp_J_kgK"]
    dT_ref   = T_hot - T0
    if dT_ref <= 0:
        raise ValueError("T_hot_K muss > initial.T_K sein (ΔT_ref > 0).")

    def _piecewise(seg_values, seg_durs):
        """Baue (times, values) als Rampe→Halten je Segment, n-fach wiederholt.
        Jedes Segment j (global) belegt exakt das Intervall
        [j·dur, (j+1)·dur] (bei konstanter Segmentdauer)."""
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
        times.append(t_now); vals.append(0.0)   # sauber auf 0 zurück
        return np.array(times), np.array(vals), t_now

    def _merge_intervals(ivs):
        """Zusammenhängende Intervalle gleicher T_inj verschmelzen."""
        merged = []
        for t0, t1, Ti in ivs:
            if merged and abs(merged[-1][1] - t0) < 1e-6 and merged[-1][2] == Ti:
                merged[-1] = (merged[-1][0], t1, Ti)
            else:
                merged.append((t0, t1, Ti))
        return merged

    # === Monatsprofil-Modus (Modus B) ===
    monthly_P = cyc.get("monthly_power_W")
    if monthly_P is not None:
        assert len(monthly_P) == 12, "cycles.monthly_power_W muss 12 Werte enthalten."
        monthly_T = cyc.get("monthly_T_inj_K")
        if monthly_T is not None:
            assert len(monthly_T) == 12, "cycles.monthly_T_inj_K muss 12 Werte enthalten."
        P = np.asarray(monthly_P, dtype=float)
        P_nom = float(np.max(np.abs(P)))
        if P_nom == 0.0:
            raise ValueError("cycles.monthly_power_W enthält nur Nullen.")
        mdot_nom  = P_nom / (cp_f * dT_ref)
        month_dur = 365.25 / 12.0 * DAY
        g_vals = (P / P_nom).tolist()                       # ∈ [−1, +1]
        # Bedarfsführung: Förder-Pumprate je Monat mit rate_mult skalieren (P<0).
        if rate_mult is not None:
            g_vals = [g * (rate_mult[m] if P[m] < 0 else 1.0)
                      for m, g in enumerate(g_vals)]
        times, vals, t_tot = _piecewise(g_vals, [month_dur] * 12)

        # Beladungs-Intervalle: Monate mit P>0 (Monat k=y·12+m belegt
        # [k·month_dur, (k+1)·month_dur]); benachbarte gleicher T_inj mergen.
        ivs = []
        for y in range(n):
            for m in range(12):
                if P[m] > 0.0:
                    k = y * 12 + m
                    T_inj = float(monthly_T[m]) if monthly_T is not None else T_hot
                    ivs.append((k * month_dur, (k + 1) * month_dur, T_inj))
        charge_intervals = _merge_intervals(ivs)

        return {
            "t_total":         t_tot,
            "cycle_power":     (times, vals),
            "charge_intervals": charge_intervals,
            "P_nom_W":         P_nom,
            "mdot_nom_kg_s":   mdot_nom,
            "mode":            "monthly",
        }

    # === 4-Phasen-Modus (Modus A) ===
    mdot_nom = cfg["operation"]["mass_flow_rate_kg_s"]
    P_nom    = mdot_nom * cp_f * dT_ref
    durs = [cyc["charge_days"] * DAY, cyc["storage_after_charge_days"] * DAY,
            cyc["discharge_days"] * DAY, cyc["storage_after_discharge_days"] * DAY]
    g_vals = [+1.0, 0.0, -1.0, 0.0]                         # Beladung/Pause/Förderung/Pause
    times, vals, t_tot = _piecewise(g_vals, durs)
    cycle_dur = sum(durs)
    charge_intervals = [(y * cycle_dur, y * cycle_dur + durs[0], T_hot)
                        for y in range(n) if durs[0] > 0]
    return {
        "t_total":         t_tot,
        "cycle_power":     (times, vals),
        "charge_intervals": charge_intervals,
        "P_nom_W":         P_nom,
        "mdot_nom_kg_s":   mdot_nom,
        "mode":            "phases",
    }


# ======================================================================
#  3) OGS-Projektdatei (axisymmetrisch)
# ======================================================================
def _se(parent, tag, text=None, **attrs):
    el = ET.SubElement(parent, tag, **{k: str(v) for k, v in attrs.items()})
    if text is not None: el.text = str(text)
    return el


def _const_prop(parent, name, value):
    p = _se(parent, "property"); _se(p, "name", name); _se(p, "type", "Constant"); _se(p, "value", value)


def _linear_T_density(parent, rho_ref, T_ref, beta):
    """Dichte ρ(T) = ρ_ref·(1 + β·(T−T_ref)) für Auftrieb (β<0)."""
    p = _se(parent, "property"); _se(p, "name", "density"); _se(p, "type", "Linear")
    _se(p, "reference_value", rho_ref)
    iv = _se(p, "independent_variable")
    _se(iv, "variable_name", "temperature")
    _se(iv, "reference_condition", T_ref)
    _se(iv, "slope", beta)


def _add_phase_fluid(phases, fluid, op):
    ph = _se(phases, "phase"); _se(ph, "type", "AqueousLiquid")
    props = _se(ph, "properties")
    beta = fluid.get("beta_1_per_K", 0.0)
    if abs(beta) > 1e-12:
        _linear_T_density(props, fluid["rho_ref_kg_m3"], fluid["T_ref_K"], beta)
    else:
        _const_prop(props, "density", fluid["rho_ref_kg_m3"])
    visc_slope = fluid.get("visc_slope_1_per_K")
    if visc_slope:
        pv = _se(props, "property"); _se(pv, "name", "viscosity"); _se(pv, "type", "Linear")
        _se(pv, "reference_value", fluid["viscosity_Pa_s"])
        iv = _se(pv, "independent_variable")
        _se(iv, "variable_name", "temperature")
        _se(iv, "reference_condition", fluid["T_ref_K"])
        _se(iv, "slope", visc_slope)
    else:
        _const_prop(props, "viscosity", fluid["viscosity_Pa_s"])
    _const_prop(props, "specific_heat_capacity", fluid["cp_J_kgK"])
    _const_prop(props, "thermal_conductivity",   fluid["lambda_W_mK"])
    _const_prop(props, "storage",                op["fluid_storage_1_per_Pa"])


def _add_phase_solid(phases, mat, op):
    ph = _se(phases, "phase"); _se(ph, "type", "Solid")
    props = _se(ph, "properties")
    _const_prop(props, "density",                mat["rho_s_kg_m3"])
    _const_prop(props, "specific_heat_capacity", mat["cp_s_J_kgK"])
    _const_prop(props, "thermal_conductivity",   mat["lambda_s_W_mK"])
    _const_prop(props, "storage",                op["solid_storage_1_per_Pa"])


def _add_medium(media, mid, mat, fluid, op, disp):
    med = _se(media, "medium", id=mid)
    phases = _se(med, "phases")
    _add_phase_fluid(phases, fluid, op)
    _add_phase_solid(phases, mat, op)
    props = _se(med, "properties")
    _const_prop(props, "porosity",     mat["porosity"])
    # Permeabilität: skalar (isotrop) oder anisotroper 2×2-Tensor (k_hor, k_ver).
    k_ver = mat.get("permeability_ver_m2")
    if k_ver is not None:
        k_hor = mat["permeability_m2"]
        _const_prop(props, "permeability", f"{k_hor:.6e} 0 0 {k_ver:.6e}")
    else:
        _const_prop(props, "permeability", mat["permeability_m2"])
    p = _se(props, "property"); _se(p, "name", "thermal_conductivity")
    _se(p, "type", "EffectiveThermalConductivityPorosityMixing")
    _const_prop(props, "thermal_longitudinal_dispersivity", disp["alpha_L_m"])
    _const_prop(props, "thermal_transversal_dispersivity",  disp["alpha_T_m"])
    _const_prop(props, "storage", 0.0)


def _curve_xml(parent, name, t, v):
    c = _se(parent, "curve"); _se(c, "name", name)
    _se(c, "coords", " ".join(f"{x:.6e}" for x in t))
    _se(c, "values", " ".join(f"{x:.6e}" for x in v))


def _const_param(parent, name, value):
    p = _se(parent, "parameter"); _se(p, "name", name); _se(p, "type", "Constant"); _se(p, "value", value)


def _curve_param(parent, name, curve_name, base):
    p = _se(parent, "parameter"); _se(p, "name", name); _se(p, "type", "CurveScaled")
    _se(p, "curve", curve_name); _se(p, "parameter", base)


def _indent(elem, level=0):
    pad = "\n" + level * "    "
    if len(elem):
        if not (elem.text and elem.text.strip()): elem.text = pad + "    "
        for child in elem: _indent(child, level + 1)
        if not (elem[-1].tail and elem[-1].tail.strip()): elem[-1].tail = pad
    if level and not (elem.tail and elem.tail.strip()): elem.tail = pad


def build_prj(cfg: dict, out_dir: Path, mesh_files: dict, curves: dict) -> Path:
    prefix = cfg["output"]["prefix"]
    fluid, op, init, sol = cfg["fluid"], cfg["operation"], cfg["initial"], cfg["solver"]
    disp = cfg["dispersion"]

    # 3D-äquivalentes Brunnenvolumen (Zylinder) für die volumetrischen Quellterme.
    # OGS integriert bei axialer Symmetrie mit 2πr, sodass
    #   ∫ q · 2πr dr dz (über die Brunnenbox) = q · V_well = Gesamtstrom.
    # ACHTUNG: Die Physical Group `hot_well_vol` ist die Brunnenspalte über die
    # VOLLE Aquiferhöhe (build_mesh: well = [S[(0, aq_j)]], und der Aquifer ist
    # dort ein einziges z-Segment). Ein Filter-Offset > 0 würde V_well
    # verkleinern, ohne das Quellgebiet zu verkleinern — der Brunnen injizierte
    # dann um den Faktor t_aq/h_screen zu viel. Deshalb: V_well aus der vollen
    # Aquiferhöhe, und Offsets werden nicht stillschweigend akzeptiert.
    _off = (cfg["well"]["screen_top_offset_m"] + cfg["well"]["screen_bottom_offset_m"])
    if abs(_off) > 1e-9:
        raise ValueError(
            "well.screen_top_offset_m / screen_bottom_offset_m werden vom 2D-Netz "
            "nicht abgebildet: die Gruppe 'hot_well_vol' umfasst immer die volle "
            "Aquiferhöhe. Ein Offset > 0 würde den Massenstrom um den Faktor "
            f"{cfg['layers']['aquifer_thickness_m']}/h_screen verfälschen. "
            "Bitte beide Offsets auf 0.0 lassen (Filter über die volle "
            "Aquiferhöhe) oder das Netz um getrennte z-Segmente erweitern."
        )
    h_screen = cfg["layers"]["aquifer_thickness_m"]
    r_well = cfg["well"]["r_well_m"]
    V_well = np.pi * r_well**2 * h_screen          # m^3 (Zylinder)

    # Basis-Amplitude der Quelle auf der p-Gleichung (cycle_power ∈[−1,1] skaliert):
    # WICHTIG: Der Quellterm der HT-Druckgleichung ist VOLUMETRISCH [m³/(m³·s)],
    # nicht massenbasiert. Der Volumenstrom ist Q = ṁ/ρ_f. Ohne Division durch ρ_f
    # würde das Wasser (und damit die Wärme) um den Faktor ρ_f (~1000) über-injiziert.
    #   q_v [1/s] -> Gesamt-Volumenstrom Q_nom bei |g|=1  (∫ q_v·2πr dV = ṁ_nom/ρ_f)
    rho_f = fluid["rho_ref_kg_m3"]
    q_mass_amp = curves["mdot_nom_kg_s"] / rho_f / V_well   # 1/s (= m³/m³/s)

    root = ET.Element("OpenGeoSysProject")

    # --- Meshes mit axially_symmetric Flag
    meshes = _se(root, "meshes")
    for key in ("domain", "top", "bottom", "far", "hot_well_vol", "hot_well_surf"):
        _se(meshes, "mesh", mesh_files[key], axially_symmetric="true")

    # --- Process
    processes = _se(root, "processes")
    proc = _se(processes, "process")
    _se(proc, "name", "HT"); _se(proc, "type", "HT"); _se(proc, "integration_order", 2)
    pv = _se(proc, "process_variables"); _se(pv, "temperature", "T"); _se(pv, "pressure", "p")
    sv = _se(proc, "secondary_variables")
    _se(sv, "secondary_variable", internal_name="darcy_velocity", output_name="darcy_velocity")
    _se(proc, "specific_body_force", "0 0" if not cfg["time"]["gravity"] else "0 -9.81")

    # --- Media
    well_mat = dict(cfg["materials"]["aquifer"])
    well_mat["permeability_m2"] = cfg["well"]["screen_permeability_m2"]
    media = _se(root, "media")
    _add_medium(media, 0, cfg["materials"]["aquifer"],        fluid, op, disp)
    _add_medium(media, 1, cfg["materials"]["caprock_top"],    fluid, op, disp)
    _add_medium(media, 2, cfg["materials"]["caprock_bottom"], fluid, op, disp)
    _add_medium(media, 3, well_mat,                           fluid, op, disp)

    # --- Time loop
    tl = _se(root, "time_loop")
    procs = _se(tl, "processes")
    p_ref = _se(procs, "process", ref="HT")
    _se(p_ref, "nonlinear_solver", "basic_picard")
    cc = _se(p_ref, "convergence_criterion")
    _se(cc, "type", "PerComponentDeltaX"); _se(cc, "norm_type", "NORM2")
    _se(cc, "reltols", f"{sol['rel_tol_T']} {sol['rel_tol_p']}")
    td = _se(p_ref, "time_discretization"); _se(td, "type", "BackwardEuler")
    ts = _se(p_ref, "time_stepping")
    _se(ts, "type", "FixedTimeStepping")
    _se(ts, "t_initial", 0.0); _se(ts, "t_end", curves["t_total"])
    steps = _se(ts, "timesteps"); pair = _se(steps, "pair")
    n_steps = int(np.ceil(curves["t_total"] / cfg["time"]["dt_seconds"]))
    _se(pair, "repeat", n_steps); _se(pair, "delta_t", cfg["time"]["dt_seconds"])

    out = _se(tl, "output")
    _se(out, "type", "VTK"); _se(out, "prefix", prefix)
    out_steps = _se(out, "timesteps"); pair = _se(out_steps, "pair")
    _se(pair, "repeat", n_steps); _se(pair, "each_steps", cfg["time"]["output_every_n_steps"])
    _se(out, "output_iteration_results", "false")
    vars_el = _se(out, "variables")
    for v in cfg["output"]["variables"]:
        _se(vars_el, "variable", v)

    # --- Parameters
    params = _se(root, "parameters")
    gradient = init.get("geothermal_gradient_K_per_m", 0.0)
    if abs(gradient) > 1e-12:
        T_surf = init.get("T_surface_K", init["T_K"])
        z_tot  = (cfg["layers"]["caprock_bottom_thickness_m"]
                  + cfg["layers"]["aquifer_thickness_m"]
                  + cfg["layers"]["caprock_top_thickness_m"])
        p_el = _se(params, "parameter")
        _se(p_el, "name", "T0")
        _se(p_el, "type", "Function")
        _se(p_el, "expression", f"{T_surf} + ({gradient:.6g})*({z_tot} - y)")
    else:
        _const_param(params, "T0", init["T_K"])
    # Druck: bei Gravitation hydrostatisch p(y)=p_top+ρ_ref·g·(z_top−y),
    # sonst konstant. z_top = Oberkante des Modells.
    if cfg["time"]["gravity"]:
        z_top = (cfg["layers"]["caprock_bottom_thickness_m"]
                 + cfg["layers"]["aquifer_thickness_m"]
                 + cfg["layers"]["caprock_top_thickness_m"])
        rho0 = fluid["rho_ref_kg_m3"]
        p_el = _se(params, "parameter")
        _se(p_el, "name", "p0"); _se(p_el, "type", "Function")
        _se(p_el, "expression", f"{init['p_Pa']:.6g} + {rho0*9.81:.6g}*({z_top:.6g} - y)")
    else:
        _const_param(params, "p0", init["p_Pa"])
    # Basis-Amplitude der Massenquelle; cycle_power ∈[−1,1] skaliert sie zeitabhängig.
    _const_param(params, "q_mass_amp", q_mass_amp)   # kg/m³/s  (ṁ_nom / V_well)
    _curve_param(params, "q_mass_well", "cycle_power", "q_mass_amp")
    # Injektionstemperatur(en): je eindeutigem Wert ein Constant-Parameter.
    tinj_values = sorted({round(iv[2], 6) for iv in curves["charge_intervals"]})
    tinj_param = {}
    for i, Ti in enumerate(tinj_values):
        name = "T_inj" if len(tinj_values) == 1 else f"T_inj_{i}"
        _const_param(params, name, Ti)
        tinj_param[Ti] = name

    # --- Curves
    cv = _se(root, "curves")
    t, v = curves["cycle_power"]; _curve_xml(cv, "cycle_power", t, v)

    # --- Process variables
    pvars = _se(root, "process_variables")

    pv_T = _se(pvars, "process_variable")
    _se(pv_T, "name", "T"); _se(pv_T, "components", 1); _se(pv_T, "order", 1)
    _se(pv_T, "initial_condition", "T0")
    bcs = _se(pv_T, "boundary_conditions")
    # Fernfeld-Thermik: T0 an Ober-/Unterkante und am Fernrand (wie MOOSE).
    for face in ("top", "bottom", "far"):
        bc = _se(bcs, "boundary_condition")
        _se(bc, "mesh", Path(mesh_files[face]).stem)
        _se(bc, "type", "Dirichlet"); _se(bc, "parameter", "T0")
    # Brunnen: Injektionstemperatur T_inj — NUR während der Beladungs-Intervalle
    # (DirichletWithinTimeInterval). Bei Förderung/Ruhe KEIN T-Zwang -> die
    # Entnahmetemperatur wird dynamisch berechnet.
    well_mesh_T = Path(mesh_files["hot_well_vol"]).stem
    for t0, t1, Ti in curves["charge_intervals"]:
        bc = _se(bcs, "boundary_condition")
        _se(bc, "mesh", well_mesh_T)
        _se(bc, "type", "DirichletWithinTimeInterval")
        _se(bc, "parameter", tinj_param[round(Ti, 6)])
        ti = _se(bc, "time_interval")
        _se(ti, "start", f"{t0:.6e}"); _se(ti, "end", f"{t1:.6e}")

    pv_p = _se(pvars, "process_variable")
    _se(pv_p, "name", "p"); _se(pv_p, "components", 1); _se(pv_p, "order", 1)
    _se(pv_p, "initial_condition", "p0")
    bcs = _se(pv_p, "boundary_conditions")
    for face in ("top", "bottom", "far"):
        bc = _se(bcs, "boundary_condition")
        _se(bc, "mesh", Path(mesh_files[face]).stem)
        _se(bc, "type", "Dirichlet"); _se(bc, "parameter", "p0")
    sts = _se(pv_p, "source_terms")
    st = _se(sts, "source_term")
    _se(st, "mesh", Path(mesh_files["hot_well_vol"]).stem)
    _se(st, "type", "Volumetric"); _se(st, "parameter", "q_mass_well")

    # --- Solvers
    nls = _se(root, "nonlinear_solvers")
    n = _se(nls, "nonlinear_solver")
    _se(n, "name", "basic_picard"); _se(n, "type", "Picard")
    _se(n, "max_iter", sol["nonlinear_iter"]); _se(n, "linear_solver", "general_linear_solver")
    lss = _se(root, "linear_solvers")
    ls = _se(lss, "linear_solver")
    _se(ls, "name", "general_linear_solver")
    eig = _se(ls, "eigen")
    # Solver konfigurierbar: "SparseLU" (direkt, robust, aber teuer bei feinem
    # Netz) oder iterativ "BiCGSTAB"/"CG" + Vorkonditionierer. Auf dem
    # strukturierten Quad-Netz konvergiert BiCGSTAB+ILUT gut und ist bei vielen
    # Unbekannten deutlich schneller.
    stype = sol.get("solver_type", "SparseLU")
    _se(eig, "solver_type", stype)
    if stype != "SparseLU":
        _se(eig, "precon_type", sol.get("precon_type", "ILUT"))
        _se(eig, "max_iteration_step", sol["linear_iter"])
        _se(eig, "error_tolerance",    sol["linear_tol"])
    _se(eig, "scaling", "true")

    _indent(root)
    prj_path = out_dir / f"{prefix}.prj"
    ET.ElementTree(root).write(prj_path, encoding="ISO-8859-1", xml_declaration=True)
    return prj_path


# ======================================================================
#  4) Run + Plots
# ======================================================================
def run_ogs(prj_path: Path) -> int:
    ogs_exe = shutil.which("ogs") or shutil.which("ogs.exe")
    if not ogs_exe:
        print("ogs.exe nicht im PATH", file=sys.stderr); return 1
    cmd = [ogs_exe, str(prj_path), "-o", str(prj_path.parent)]
    print(">>", " ".join(cmd))
    return subprocess.call(cmd)


def make_plots(cfg: dict, out_dir: Path, rate_mult=None) -> None:
    try:
        import pyvista as pv
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("pyvista/matplotlib fehlt, Plots uebersprungen."); return

    # Plots in den konventionellen figures/-Ordner der Übung (nicht out/).
    # Neben den LAUF, nicht neben das Skript: sonst ueberschreiben sich zwei
    # Laeufe gegenseitig die Bilder und sie liegen nicht bei ihren Daten.
    figdir = Path(out_dir) / "figures"; figdir.mkdir(parents=True, exist_ok=True)
    prefix = cfg["output"]["prefix"]
    files = sorted(out_dir.glob(f"{prefix}_ts_*_t_*.vtu"),
                   key=lambda p: int(re.search(r"_ts_(\d+)_", p.name).group(1)))
    if not files: print("Keine VTU-Ausgaben."); return

    T0    = cfg["initial"]["T_K"]
    T_hot = cfg["operation"]["T_hot_K"]
    cp_f  = cfg["fluid"]["cp_J_kgK"]
    z_base = cfg["domain"]["z_base_m"]
    t_cb   = cfg["layers"]["caprock_bottom_thickness_m"]
    t_aq   = cfg["layers"]["aquifer_thickness_m"]
    z_aq_bot = z_base + t_cb
    z_aq_mid = z_aq_bot + t_aq / 2.0
    z_aq_top = z_aq_bot + t_aq

    monthly_P = cfg["cycles"].get("monthly_power_W")
    monthly   = monthly_P is not None
    n_cyc     = cfg["cycles"]["n_cycles"]
    YEAR      = 365.25
    month_dur = YEAR / 12.0
    MONTHS = ["Jan","Feb","Mär","Apr","Mai","Jun","Jul","Aug","Sep","Okt","Nov","Dez"]

    def day_of(f):
        return int(re.search(r"_t_(\d+)", f.name).group(1)) / DAY

    # --- rcp-Feld (axisymm. Volumen) für Energiebilanz ---
    domain = pv.read(out_dir / f"{prefix}_domain.vtu")
    mid = domain.cell_data["MaterialIDs"]
    cell_area = domain.compute_cell_sizes()["Area"]
    r_centroid = domain.cell_centers().points[:, 0]
    cell_vol_axi = 2 * np.pi * r_centroid * cell_area
    a = cfg["materials"]["aquifer"]; c = cfg["materials"]["caprock_top"]
    rcp_aq = a["porosity"]*1000*cp_f + (1-a["porosity"])*a["rho_s_kg_m3"]*a["cp_s_J_kgK"]
    rcp_cr = c["porosity"]*1000*cp_f + (1-c["porosity"])*c["rho_s_kg_m3"]*c["cp_s_J_kgK"]
    rcp = np.where(mid == 0, rcp_aq, np.where(np.isin(mid, [1, 2]), rcp_cr, rcp_aq))

    # --- Beobachtungspunkte (Front ~ einige 10er m bei MOOSE-Skala) ---
    probes = {
        "Brunnen (r≈0.3 m)": (0.3,  z_aq_mid, 0),
        "r=5 m":             (5.0,  z_aq_mid, 0),
        "r=15 m":            (15.0, z_aq_mid, 0),
        "r=30 m":            (30.0, z_aq_mid, 0),
        "Cap Rock (+5 m)":   (2.0,  z_aq_top + 5.0, 0),
    }
    probe_pts = pv.PolyData(np.array(list(probes.values())))

    # === EINZELDURCHLAUF: Sonde + Energie je Datei ===
    times = []
    series = {k: [] for k in probes}
    E_tot = []
    for f in files:
        m = pv.read(f)
        times.append(day_of(f))
        s = probe_pts.sample(m)
        for i, k in enumerate(probes):
            series[k].append(float(s["T"][i]))
        cT = m.point_data_to_cell_data(pass_point_data=False)["T"]
        E_tot.append(float(np.sum(rcp * cell_vol_axi * (cT - T0))) / 1e9)
    times = np.array(times); E_tot = np.array(E_tot)

    # === Nominal-Energiebilanz (Soll) aus Monatsprofil ===
    def power_at(day):
        if not monthly:
            return 0.0
        m = int((day % YEAR) // month_dur) % 12
        return float(monthly_P[m])
    P_series = np.array([power_at(d) for d in times])
    # kumulierte Soll-Energie ∫P dt
    E_nom = np.zeros_like(times)
    for i in range(1, len(times)):
        E_nom[i] = E_nom[i-1] + 0.5*(P_series[i]+P_series[i-1])*(times[i]-times[i-1])*DAY/1e9

    # ---------------------------------------------------------------
    #  1) T-Feld-Montage: repräsentatives (spätes) Jahr, 6 Monate
    # ---------------------------------------------------------------
    Y = max(0, n_cyc - 5) if monthly else 0        # eingeschwungenes Jahr
    snap_months = [1, 3, 5, 7, 9, 11]              # Feb, Apr, Jun, Aug, Okt, Dez
    keypoints = [((Y*YEAR + (mo+1)*month_dur), f"{MONTHS[mo]} J{Y+1}")
                 for mo in snap_months]
    r_view = 60.0; z_lo = z_aq_bot - 22.0; z_hi = z_aq_top + 22.0
    plotter = pv.Plotter(off_screen=True, window_size=(1800, 640),
                         shape=(1, len(keypoints)), border=True)
    pv.OFF_SCREEN = True
    for i, (day, label) in enumerate(keypoints):
        f = min(files, key=lambda p: abs(day_of(p) - day))
        m = pv.read(f).clip_box((0, r_view, z_lo, z_hi, -1, 1), invert=False)
        # zur vollen (r,z)-Sektion spiegeln (axialsymmetrisch)
        m = m.merge(m.reflect((1, 0, 0), point=(0, 0, 0)))
        plotter.subplot(0, i)
        plotter.add_mesh(m, scalars="T", cmap="coolwarm", clim=[T0, T_hot],
                         show_scalar_bar=(i == len(keypoints)-1),
                         scalar_bar_args={"title": "T [K]", "vertical": True,
                                          "position_x": 0.05} if i == len(keypoints)-1 else None,
                         show_edges=False)
        plotter.add_text(label, font_size=11, position="upper_edge")
        plotter.view_xy(); plotter.camera.zoom(1.5)
    plotter.screenshot(str(figdir / "T_field_2D.png"))
    plotter.close()
    print("  saved T_field_2D.png")

    # ---------------------------------------------------------------
    #  2) T(t) an Beobachtungspunkten (+ Monatsleistung als Overlay)
    #     2 Reihen: gesamt + Zoom auf die letzten 3 Jahre
    # ---------------------------------------------------------------
    def plot_T(ax, tmin, tmax, show_leg):
        axp = ax.twinx()
        if monthly:
            axp.step(times, P_series/1e3, where="post", color="0.6", lw=0.8, alpha=0.6)
            axp.axhline(0, color="0.6", lw=0.5, ls=":")
            axp.set_ylabel("Monatsleistung P [kW]", color="0.5", fontsize=9)
            axp.tick_params(axis="y", labelcolor="0.5")
        for k, T in series.items():
            ax.plot(times, np.array(T), lw=1.3, label=k, zorder=3)
        ax.axhline(T0,    color="k", lw=0.6, ls=":")
        ax.axhline(T_hot, color="k", lw=0.6, ls=":")
        ax.set_xlim(tmin, tmax); ax.set_ylim(T0-6, T_hot+6)
        ax.set_xlabel("Zeit [Tage]"); ax.set_ylabel("T [K]")
        ax.set_zorder(axp.get_zorder()+1); ax.patch.set_visible(False)
        ax.grid(alpha=0.3)
        if show_leg: ax.legend(loc="upper left", fontsize=8, ncol=3)

    fig, axs = plt.subplots(2, 1, figsize=(13, 8))
    plot_T(axs[0], 0, times[-1], True)
    axs[0].set_title("ATES 2D radial — T(t) an Beobachtungspunkten (Brunnentemp. = dynamische Entnahmetemperatur)")
    ny_zoom = min(3, n_cyc)
    z0 = max(0, times[-1] - ny_zoom*YEAR)
    plot_T(axs[1], z0, times[-1], False)
    axs[1].set_title("Zoom: gesamter Zeitraum" if ny_zoom >= n_cyc
                     else f"Zoom: letzte {ny_zoom} Betriebsjahre")
    fig.tight_layout()
    fig.savefig(figdir / "T_vs_time.png", dpi=130); plt.close(fig)
    print("  saved T_vs_time.png")

    # ---------------------------------------------------------------
    #  2b) Vollständige Brunnentemperatur-Kurve über den Zeitraum
    #      (= Injektionstemperatur bei Beladung / dynamische Entnahme-
    #       temperatur bei Förderung). Daten auch in _timeseries.csv.
    # ---------------------------------------------------------------
    T_inj = cfg["operation"]["T_hot_K"]
    well = np.array(series[list(probes.keys())[0]])
    fig, ax = plt.subplots(figsize=(13, 4.6))
    if monthly:
        for y in range(n_cyc):
            for mo in range(12):
                t0m = y * YEAR + mo * month_dur
                if monthly_P[mo] > 0:
                    ax.axvspan(t0m, t0m + month_dur, color="red",  alpha=0.06)
                elif monthly_P[mo] < 0:
                    ax.axvspan(t0m, t0m + month_dur, color="blue", alpha=0.06)
    ax.plot(times, well - 273.15, color="#b22222", lw=1.6, zorder=3,
            label="Brunnentemperatur")
    ax.axhline(T_inj - 273.15, color="k", lw=0.8, ls="--", label=f"T_inj = {T_inj-273.15:.0f} °C")
    ax.axhline(T0 - 273.15,   color="k", lw=0.8, ls=":",  label=f"T_amb = {T0-273.15:.0f} °C")
    ax.set_xlabel("Zeit [Tage]"); ax.set_ylabel("Brunnentemperatur [°C]")
    ax.set_title("ATES 2D radial — Brunnentemperatur über den gesamten Zeitraum\n"
                 "(Beladung/rot: Injektion bei T_inj — Förderung/blau: dynamische Entnahmetemperatur)")
    ax.grid(alpha=0.3); ax.legend(loc="best", fontsize=9)
    ax.set_xlim(0, times[-1]); ax.set_ylim(T0 - 273.15 - 4, T_inj - 273.15 + 4)
    fig.tight_layout()
    fig.savefig(figdir / "well_temperature.png", dpi=130); plt.close(fig)
    print("  saved well_temperature.png")

    # ---------------------------------------------------------------
    #  3) Energiebilanz: gespeicherte Energie + Soll (kumuliert)
    # ---------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(13, 4.8))
    # Beladungs-/Fördermonate als leichte Hintergrundschattierung (nur 1. Reihe Jahre sparsam)
    if monthly:
        for y in range(n_cyc):
            for mo in range(12):
                if monthly_P[mo] > 0:
                    t0m = y*YEAR + mo*month_dur
                    ax.axvspan(t0m, t0m+month_dur, color="red",  alpha=0.05)
                elif monthly_P[mo] < 0:
                    t0m = y*YEAR + mo*month_dur
                    ax.axvspan(t0m, t0m+month_dur, color="blue", alpha=0.05)
    ax.plot(times, E_tot, "k-", lw=1.8, label="Gespeicherte Wärme (Simulation)")
    ax.set_xlabel("Zeit [Tage]"); ax.set_ylabel("Wärme über T0 [GJ]")
    ax.set_title("ATES 2D radial — gespeicherte Wärme im Untergrund")
    ax.grid(alpha=0.3); ax.set_xlim(0, times[-1]); ax.set_ylim(bottom=0)
    # Jahres-Recovery (letztes volles Jahr): (E_max−E_min)/E_max im Jahr
    if monthly and n_cyc >= 1:
        yl0, yl1 = (n_cyc-1)*YEAR, n_cyc*YEAR
        seg = (times >= yl0) & (times <= yl1)
        if seg.sum() > 2:
            Emax, Emin = E_tot[seg].max(), E_tot[seg].min()
            if Emax > 0:
                rec = (Emax - Emin)/Emax*100
                ax.annotate(f"Jahres-Entladung (letztes Jahr): {rec:.0f} %\n"
                            f"ΔE = {Emax-Emin:.0f} GJ",
                            xy=(0.98, 0.05), xycoords="axes fraction",
                            ha="right", va="bottom", fontsize=10,
                            bbox=dict(boxstyle="round", fc="w", alpha=0.8))
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(figdir / "energy_balance.png", dpi=130); plt.close(fig)
    print("  saved energy_balance.png")

    # ---------------------------------------------------------------
    #  4) Recovery-Effizienz R pro Jahr (MOOSE-Kernmetrik)
    #     R = (T̄_prod − T_amb)/(T_inj − T_amb),  T̄_prod = |P|-gewichtetes
    #     Mittel der Brunnentemperatur über die Fördermonate.
    # ---------------------------------------------------------------
    if monthly:
        T_inj = cfg["operation"]["T_hot_K"]; T_amb = T0
        well = np.array(series["Brunnen (r≈0.3 m)"])
        R_years, Tp_years, yy = [], [], []
        for y in range(n_cyc):
            num = den = 0.0
            for mo in range(12):
                if monthly_P[mo] < 0:
                    t0m, t1m = y*YEAR + mo*month_dur, y*YEAR + (mo+1)*month_dur
                    seg = (times >= t0m) & (times < t1m)
                    if seg.any():
                        w = abs(monthly_P[mo])
                        num += w * well[seg].mean(); den += w
            if den > 0:
                Tp = num/den
                Tp_years.append(Tp); R_years.append((Tp - T_amb)/(T_inj - T_amb)); yy.append(y+1)
        if R_years:
            fig, ax = plt.subplots(figsize=(10, 4.4))
            ax.bar(yy, np.array(R_years)*100, color="#c0504d", alpha=0.85, width=0.7)
            for x, R in zip(yy, R_years):
                ax.text(x, R*100+1, f"{R*100:.0f}%", ha="center", va="bottom", fontsize=8)
            ax.set_xlabel("Betriebsjahr"); ax.set_ylabel("Temperaturhub-Ausnutzung [%]")
            ax.set_title("ATES 2D radial — Ausnutzung des Temperaturhubs\n"
                         f"R = (T̄_Förder − T_amb)/(T_inj − T_amb),  T_inj={T_inj-273.15:.0f} °C, "
                         f"T_amb={T_amb-273.15:.0f} °C")
            ax.set_ylim(0, max(100, max(R_years)*100*1.15)); ax.grid(alpha=0.3, axis="y")
            # Ohne diesen Hinweis stehen zwei verschieden definierte
            # "Rueckgewinnungsgrade" fuer denselben Lauf nebeneinander:
            # hier temperaturbasiert (~16 %), im Pruefblatt energetisch
            # (~53 %). Das hat schon einmal fuer Ratlosigkeit gesorgt.
            ax.text(0.5, -0.30, "Das ist NICHT der energetische Rueckgewinnungsgrad eta = E_aus/E_ein "
                    "- der steht im Pruefblatt (0_pruefblatt.png) und liegt deutlich hoeher.",
                    transform=ax.transAxes, ha="center", fontsize=9, color="#5a5a5a")
            ax.axhline(100, color="k", lw=0.5, ls=":")
            fig.tight_layout()
            fig.savefig(figdir / "temperaturhub_ausnutzung.png", dpi=130); plt.close(fig)
            print("  saved temperaturhub_ausnutzung.png")
            print("  Temperaturhub-Ausnutzung pro Jahr [%]:", [f"{r*100:.0f}" for r in R_years])

        # -----------------------------------------------------------
        #  5) Bedarfsführung: gelieferte vs. geforderte Förderleistung
        # -----------------------------------------------------------
        if rate_mult is not None:
            mult = rate_mult
            dT_ref = T_inj - T_amb
            mo_of = (np.floor((times % YEAR) / month_dur).astype(int)) % 12
            P_arr = np.array([monthly_P[m] for m in mo_of], dtype=float)
            fac   = np.array([mult[m] for m in mo_of])
            is_prod = P_arr < 0
            # gelieferte Entnahmeleistung = ṁ·c_p·(T_prod−T_amb) = |P|·fac·(T_w−T_amb)/ΔT_ref
            deliver = np.where(is_prod,
                               np.abs(P_arr) * fac * (well - T_amb) / dT_ref, 0.0)
            demand_p = np.where(is_prod, np.abs(P_arr), 0.0)
            fig, ax = plt.subplots(figsize=(13, 4.6))
            ny_d = min(4, n_cyc)
            z0 = max(0, times[-1] - ny_d*YEAR)
            sel = times >= z0
            ax.plot(times[sel], demand_p[sel]/1e3, color="0.4", lw=1.6, ls="--",
                    label="Geforderte Förderleistung |P|")
            ax.plot(times[sel], deliver[sel]/1e3, color="#1f77b4", lw=1.6,
                    label="Gelieferte Förderleistung (bedarfsgeführt)")
            ax.set_xlabel("Zeit [Tage]"); ax.set_ylabel("Förderleistung [kW]")
            span = "gesamter Zeitraum" if ny_d >= n_cyc else f"letzte {ny_d} Jahre"
            ax.set_title("ATES 2D radial — Bedarfsgeführte Förderung: gelieferte vs. geforderte Leistung "
                         f"({span})")
            ax.grid(alpha=0.3); ax.legend(loc="upper right", fontsize=9)
            ax.set_ylim(bottom=0)
            fig.tight_layout()
            fig.savefig(figdir / "demand_power.png", dpi=130); plt.close(fig)
            print("  saved demand_power.png")
            print("  konvergierte Ratenfaktoren:", [round(x, 2) for x in mult])

    # ---------------------------------------------------------------
    #  CSV-Ausgabe (Zeitreihe + Jahres-Kennzahlen), analog MOOSE
    # ---------------------------------------------------------------
    import csv
    keys = list(probes.keys())
    # Förder-/Bedarfsleistung [kW] (auch bei fester Rate mit mult=1)
    if monthly:
        T_inj = cfg["operation"]["T_hot_K"]
        mult = rate_mult if rate_mult is not None else [1.0] * 12
        mo_of = (np.floor((times % YEAR) / month_dur).astype(int)) % 12
        P_arr = np.array([monthly_P[m] for m in mo_of], float)
        fac = np.array([mult[m] for m in mo_of])
        prod = P_arr < 0
        well_a = np.array(series[keys[0]])
        demand_kw  = np.where(prod, np.abs(P_arr), 0.0) / 1e3
        deliver_kw = np.where(prod, np.abs(P_arr) * fac * (well_a - T0) / (T_inj - T0), 0.0) / 1e3

    ts_csv = out_dir / f"{prefix}_timeseries.csv"
    with open(ts_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        hdr = ["day"] + [f"T_{k}_K" for k in ["well", "r5", "r15", "r30", "caprock"]] + ["E_stored_GJ"]
        if monthly:
            hdr += ["P_demand_kW", "P_delivered_kW"]
        w.writerow(hdr)
        for i in range(len(times)):
            row = [f"{times[i]:.3f}"] + [f"{series[k][i]:.4f}" for k in keys] + [f"{E_tot[i]:.4f}"]
            if monthly:
                row += [f"{demand_kw[i]:.3f}", f"{deliver_kw[i]:.3f}"]
            w.writerow(row)
    print(f"  saved {ts_csv.name}")

    if monthly:
        sm_csv = out_dir / f"{prefix}_summary.csv"
        with open(sm_csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["year", "recovery_R_percent", "T_prod_mean_degC"])
            for yr, R, Tp in zip(yy, R_years, Tp_years):
                w.writerow([yr, f"{R*100:.1f}", f"{Tp-273.15:.2f}"])
        print(f"  saved {sm_csv.name}")
        # Konvergierte monatliche Förder-Ratenfaktoren (Bedarfsführung)
        if rate_mult is not None:
            rc_csv = out_dir / f"{prefix}_prod_rate_factors.csv"
            with open(rc_csv, "w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh)
                w.writerow(["month_index", "month", "P_W", "rate_factor"])
                for mo in range(12):
                    if monthly_P[mo] < 0:
                        w.writerow([mo, MONTHS[mo], monthly_P[mo], f"{rate_mult[mo]:.3f}"])
            print(f"  saved {rc_csv.name}")


# ======================================================================
#  Bedarfsgeführte Förderung: Fördertemperatur messen
# ======================================================================
def measure_produced_T(cfg: dict, out_dir: Path) -> dict:
    """|P|-gewichtete mittlere Brunnentemperatur je Fördermonat (P<0),
    gemittelt über die eingeschwungenen Jahre (zweite Laufhälfte)."""
    import pyvista as pv
    prefix = cfg["output"]["prefix"]
    files = sorted(out_dir.glob(f"{prefix}_ts_*_t_*.vtu"),
                   key=lambda p: int(re.search(r"_ts_(\d+)_", p.name).group(1)))
    if not files:
        return {}
    z_aq_mid = (cfg["domain"]["z_base_m"]
                + cfg["layers"]["caprock_bottom_thickness_m"]
                + cfg["layers"]["aquifer_thickness_m"] / 2.0)
    well_pt = pv.PolyData(np.array([[0.3, z_aq_mid, 0.0]]))
    YEAR = 365.25; month_dur = YEAR / 12.0
    n_cyc = cfg["cycles"]["n_cycles"]
    monthly_P = cfg["cycles"]["monthly_power_W"]
    y_start = max(1, n_cyc // 2)                 # eingeschwungene Jahre
    days, Twell = [], []
    for f in files:
        d = int(re.search(r"_t_(\d+)", f.name).group(1)) / DAY
        if d < y_start * YEAR:
            continue
        m = pv.read(f)
        Twell.append(float(well_pt.sample(m)["T"][0])); days.append(d)
    days = np.array(days); Twell = np.array(Twell)
    result = {}
    for mo in range(12):
        if monthly_P[mo] >= 0:
            continue
        vals = []
        for y in range(y_start, n_cyc):
            t0m, t1m = y * YEAR + mo * month_dur, y * YEAR + (mo + 1) * month_dur
            seg = (days >= t0m) & (days < t1m)
            if seg.any():
                vals.append(Twell[seg].mean())
        if vals:
            result[mo] = float(np.mean(vals))
    return result


# ======================================================================
#  Hauptfunktion
# ======================================================================
def main() -> int:
    ap = argparse.ArgumentParser(description="ATES 2D radial demo")
    ap.add_argument("--no-mesh",  action="store_true")
    ap.add_argument("--no-run",   action="store_true")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--years", type=int, default=None,
                    help="Override cycles.n_cycles (Betriebsjahre) für schnelle Tests.")
    ap.add_argument("--dt-days", type=float, default=None,
                    help="Override time.dt_seconds (in Tagen).")
    ap.add_argument("--production-control", choices=["demand", "fixed"], default=None,
                    help="Override well.production_control.")
    args = ap.parse_args()

    if args.years is not None:
        CONFIG["cycles"]["n_cycles"] = args.years
    if args.dt_days is not None:
        CONFIG["time"]["dt_seconds"] = args.dt_days * DAY
    if args.production_control is not None:
        CONFIG["well"]["production_control"] = args.production_control

    out_dir = Path(CONFIG["output"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = CONFIG["output"]["prefix"]
    msh_path = out_dir / f"{prefix}.msh"

    if not args.no_mesh:
        print("[1/4] gmsh: 2D Radialmesh erzeugen ...")
        build_mesh(CONFIG, out_dir)
        print(f"      -> {msh_path}")
        print("[2/4] msh2vtu ...")
        mesh_files = convert_mesh(CONFIG, msh_path, out_dir)
    else:
        mesh_files = {k: f"{prefix}_physical_group_{n}.vtu" if "physical" in n else f"{prefix}_{n}.vtu"
                      for k, n in [("domain","domain"),("top","top"),("bottom","bottom"),
                                   ("far","far"),("hot_well_vol","hot_well_vol"),
                                   ("hot_well_surf","hot_well_surf")]}
        # Fix prefix
        mesh_files = {
            "domain":         f"{prefix}_domain.vtu",
            "top":            f"{prefix}_physical_group_top.vtu",
            "bottom":         f"{prefix}_physical_group_bottom.vtu",
            "far":            f"{prefix}_physical_group_far.vtu",
            "hot_well_vol":   f"{prefix}_physical_group_hot_well_vol.vtu",
            "hot_well_surf": f"{prefix}_physical_group_hot_well_surf.vtu",
        }

    # --- Förderregelung: fest oder bedarfsgeführt (iterativ) ---
    monthly = CONFIG["cycles"].get("monthly_power_W") is not None
    demand  = (CONFIG["well"].get("production_control", "fixed") == "demand"
               and monthly and not args.no_run)
    T_amb  = CONFIG["initial"]["T_K"]
    dT_ref = CONFIG["operation"]["T_hot_K"] - T_amb
    max_fac = CONFIG["well"].get("max_rate_factor", 6.0)
    monthly_P = CONFIG["cycles"].get("monthly_power_W")

    rate_mult = [1.0] * 12
    n_iter = CONFIG["well"].get("demand_iterations", 3) if demand else 1
    for it in range(n_iter):
        tag = f" (Bedarfs-Iteration {it+1}/{n_iter})" if demand else ""
        print(f"[3/4] OGS-Projektdatei ...{tag}")
        curves = build_cycle_curves(CONFIG, rate_mult=rate_mult if demand else None)
        prj_path = build_prj(CONFIG, out_dir, mesh_files, curves)
        print(f"      -> {prj_path}  (t_end = {curves['t_total']/DAY:.0f} d)")
        if args.no_run:
            return 0
        print("[4/4] OGS starten ...")
        rc = run_ogs(prj_path)
        if rc != 0:
            return rc
        if demand and it < n_iter - 1:
            Tp = measure_produced_T(CONFIG, out_dir)
            new_mult = list(rate_mult)
            for mo, T in Tp.items():
                target = min(dT_ref / max(T - T_amb, 1.0), max_fac)   # halte Leistung
                new_mult[mo] = 0.5 * rate_mult[mo] + 0.5 * target      # Unterrelaxation
            print("  T_prod_avg [degC]:", {mo: round(T-273.15, 1) for mo, T in Tp.items()})
            print("  rate_factors:", [round(x, 2) for x in new_mult])
            rate_mult = new_mult

    if not args.no_plots:
        print("[5/4] Plots erzeugen ...")
        make_plots(CONFIG, out_dir, rate_mult=(rate_mult if demand else None))

    # --- automatischer Pruef- und Auswertebericht (ates_report.py) ----
    # Das Modul liegt in exercises/ates/. Schlaegt der Bericht fehl, darf der
    # Lauf davon nicht betroffen sein -> auto_report faengt selbst alles ab.
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

# --- Temperaturen -----------------------------------------------------
# Die Aquifertemperatur ist zugleich die Referenztemperatur des Fluids.
# Das ist kein Zufall, sondern Absicht: OGS rechnet die Viskositaet als
# mu(T) = mu_ref * (1 + slope*(T - T_ref)). Mit T_ref = T_aquifer gilt am
# ungestoerten Aquifer exakt mu = fluid.viscosity_Pa_s, und die Umrechnung
# kf -> k weiter unten ist ohne Umweg richtig.
T_AQ_K = FALL["T_aquifer_C"] + 273.15
C["initial"]["T_K"] = T_AQ_K
C["fluid"]["T_ref_K"] = T_AQ_K
C["operation"]["T_hot_K"] = FALL["T_injektion_C"] + 273.15
# operation.T_cold_K wird hier ABSICHTLICH nicht gesetzt: kein Modell liest
# den Schluessel. Bei Foerderung gibt es ueberhaupt keine Temperatur-
# Randbedingung (nur DirichletWithinTimeInterval waehrend der Beladung), die
# Foerdertemperatur ist rein dynamisch. Ihn zu setzen wuerde eine Steuerung
# vortaeuschen, die es nicht gibt.

MU = C["fluid"]["viscosity_Pa_s"]        # bei T_ref, hier also bei T_aquifer

# --- Plausibilitaet der Eingaben --------------------------------------
# OGS rechnet die Viskositaet als GERADE: mu(T) = mu_ref*(1 + slope*(T-T_ref)).
# Eine Gerade hat eine Nullstelle. Bei slope = -1.28e-2 liegt sie
# 1/0.0128 = 78.1 K ueber der Aquifertemperatur, hier also bei 88.1 GradC.
# Darueber wird die Viskositaet NEGATIV und das Modell rechnet stillschweigend
# Unsinn - es gibt dafuer keine Pruefung in OGS und keine Zeile im Pruefblatt.
# Schon ab etwa 80 GradC ist die Gerade deutlich zu flach. Wer HT-ATES rechnen
# will, braucht ein anderes Viskositaetsmodell, nicht nur ein hoeheres T_inj.
_T_NULL = FALL["T_aquifer_C"] - 1.0 / C["fluid"]["visc_slope_1_per_K"]
if FALL["T_injektion_C"] > _T_NULL - 10.0:
    raise ValueError(
        f"T_injektion_C = {FALL['T_injektion_C']:.0f} GradC ist zu hoch fuer "
        f"das lineare Viskositaetsmodell: mu(T) wird bei "
        f"{_T_NULL:.1f} GradC null und darueber negativ.\n"
        f"Zulaessig ist hier bis {_T_NULL - 10.0:.1f} GradC. Fuer hoehere "
        f"Injektionstemperaturen muss fluid.visc_slope_1_per_K angepasst "
        f"werden (oder auf None gesetzt, dann rechnet OGS mit konstanter "
        f"Viskositaet - was den Auftrieb unterschaetzt).")


def _k(kf_m_s: float) -> float:
    """Durchlaessigkeitsbeiwert kf [m/s] -> Permeabilitaet k [m2].

        k = kf * mu / (rho_f * g)

    Das mu ist NICHT frei waehlbar - es muss dasjenige sein, mit dem OGS
    bei Aquifertemperatur rechnet, sonst hat das Modell eine andere
    Durchlaessigkeit als der Standort. Genau das war in der
    Vorgaengerfassung passiert: dort stand k = 6.12e-11 m2 mit dem
    Kommentar "kf = 6e-4 m/s", umgerechnet war aber mit mu = 1.0e-3
    (Wasser bei 20 GradC) worden, waehrend das Modell bei 10 GradC mit
    mu = 1.3e-3 rechnet. Das Modell war damit 30 % zu dicht:
    K = 6.12e-11 * 9810 / 1.3e-3 = 4.62e-4 m/s statt 6.0e-4 m/s.
    """
    return kf_m_s * MU / (RHO_F * G)


# --- Aquifer ----------------------------------------------------------
_a = FALL["aquifer"]
C["layers"]["aquifer_thickness_m"] = _a["maechtigkeit_m"]
aq = C["materials"]["aquifer"]
aq["permeability_m2"] = _k(_a["kf_m_s"])
aq["permeability_ver_m2"] = _k(_a["kf_m_s"])          # isotrop
aq["porosity"] = _a["porositaet"]
aq["rho_s_kg_m3"] = _a["dichte_korn_kg_m3"]
aq["cp_s_J_kgK"] = _a["waermekapazitaet_korn_J_kgK"]
aq["lambda_s_W_mK"] = _a["waermeleitfaehigkeit_korn_W_mK"]

# --- Deckgestein ------------------------------------------------------
_d = FALL["deckgestein"]
C["layers"]["caprock_top_thickness_m"] = _d["maechtigkeit_m"]
C["layers"]["caprock_bottom_thickness_m"] = _d["maechtigkeit_m"]
for _key in ("caprock_top", "caprock_bottom"):
    cr = C["materials"][_key]
    cr["permeability_m2"] = _d["permeabilitaet_m2"]
    cr["permeability_ver_m2"] = _d["permeabilitaet_m2"] * _d["anisotropie_vertikal"]
    cr["porosity"] = _d["porositaet"]
    cr["rho_s_kg_m3"] = _d["dichte_korn_kg_m3"]
    cr["cp_s_J_kgK"] = _d["waermekapazitaet_korn_J_kgK"]
    cr["lambda_s_W_mK"] = _d["waermeleitfaehigkeit_korn_W_mK"]

# --- Brunnen ----------------------------------------------------------
C["well"]["r_well_m"] = FALL["brunnenradius_m"]
C["well"]["screen_top_offset_m"] = 0.0        # Filter ueber die volle Hoehe
C["well"]["screen_bottom_offset_m"] = 0.0     # (2D-Netz kann es nicht anders)
C["well"]["screen_permeability_m2"] = _k(_a["kf_m_s"])
# "fixed": der Massenstrom folgt dem Lastprofil. Die Vorgabe des Modells
# waere "demand" - dabei regelt der Brunnen die Foerderrate nach, bis der
# Bedarf gedeckt ist. Das ist hier ausdruecklich NICHT gewollt: eine
# nachgeregelte Rate verhindert, dass sich der Untergrund auflaedt, und
# druckt den Rueckgewinnungsgrad von 53 % auf 13 %.
C["well"]["production_control"] = "fixed"

# --- Lastprofil und Massenstrom --------------------------------------
P = FALL["monatsleistung_W"]
C["cycles"]["monthly_power_W"] = P
C["cycles"]["monthly_T_inj_K"] = None
C["cycles"]["n_cycles"] = FALL["betriebsjahre"]
C["cycles"]["ramp_days"] = 3.0
# Nur der Vollstaendigkeit halber gesetzt: im Monatsprofil-Modus liest das
# Modell diesen Schluessel NICHT, es rechnet mdot = P_max/(c_f*dT_ref) mit
# dT_ref = T_hot_K - initial.T_K. Hier wird genau diese Formel eingetragen,
# damit die Konsolenausgabe nicht etwas anderes behauptet als die Rechnung.
DT_REF = FALL["T_injektion_C"] - FALL["T_aquifer_C"]
C["operation"]["mass_flow_rate_kg_s"] = (
    max(abs(p) for p in P) / (C["fluid"]["cp_J_kgK"] * DT_REF))

# --- Physik, die nicht zum Standort gehoert ---------------------------
C["domain"]["r_max_m"] = 1000.0
C["time"]["gravity"] = True          # Auftrieb an. Bei dieser Durchlaessigkeit
                                     # ist er der Hauptverlustpfad, nicht die
                                     # Waermeleitung: Ra = 530 gegen kritisch 40.
C["dispersion"]["alpha_L_m"] = 2.0   # identisch zu Fall 2, damit die beiden
C["dispersion"]["alpha_T_m"] = 0.2   # Faelle vergleichbar bleiben

# Thermische Ausdehnung: rho(T) = rho_ref * (1 + beta*(T - T_ref)). Das ist
# die EINZIGE Zahl, die den Auftrieb steuert - und der Auftrieb ist hier der
# Hauptverlustpfad (Ra ~ 440 gegen kritisch 40), er begruendet die 120 m
# Deckgestein. Explizit gesetzt, damit er nicht stillschweigend geerbt wird.
# Es ist eine Linearisierung um T_ref: rho(60 GradC) = 980 statt real 983,2 -
# die treibende Dichtedifferenz ist rund 20 % zu gross. Qualitativ aendert das
# nichts, quantitativ ist der Auftriebsverlust systematisch etwas ueberschaetzt.
C["fluid"]["beta_1_per_K"] = -4.0e-4

# Iterationsbudget der Picard-Schleife. Die Vorgabe des Modells ist 20 - und
# genau hier lauert eine Falle: das 2D-Modell schreibt FixedTimeStepping, kann
# den Zeitschritt bei Nichtkonvergenz also NICHT verkleinern (das 3D-Modell
# kann es). Gleichzeitig rechnet dieser Fall 30 statt 2 Betriebsjahre. Der
# ausgelieferte Parametersatz braucht hoechstens 8 Iterationen, aber jede
# Verschaerfung (hoehere Injektionstemperatur, groesseres kf, groessere
# Monatsleistungen) kann das kippen - und der Lauf stirbt dann nach Stunden.
C["solver"]["nonlinear_iter"] = 50

# Speicherkoeffizient. OGS-HT liest `storage` AUSSCHLIESSLICH aus der
# Solid-Phase - der Fluid-Wert im CONFIG ist tote Konfiguration. Hier
# explizit gesetzt, damit er nicht stillschweigend aus dem Uebungsskript
# geerbt wird:  storage = S_s / (rho g), S_s = 1e-5 1/m (gespannter
# Aquifer, Lehrbuch) -> 1.0e-9 1/Pa.
# Die Ergebnisse haengen praktisch nicht daran: die hydraulische
# Diffusivitaet K/S_s liegt zwischen 4.6 und 460 m2/s, ein Drucksignal
# durchquert 150 m also in Sekunden bis Minuten - immer weit unter dem
# Zeitschritt von einem Tag. Das Druckfeld ist quasistationaer. Belegt an
# zwei fertigen Laeufen mit Faktor 100 Unterschied im Speicher-
# koeffizienten: beide treffen die stationaere Thiem-Vorhersage auf 6 %.
C["operation"]["solid_storage_1_per_Pa"] = 1.0e-9

# --- Netz -------------------------------------------------------------
# Vierteilig radial: Brunnen | Nahfeld | Fahnenband | Fernfeld. Mit EINER
# geometrischen Progression vom Brunnen bis zum Fahnenrand wird entweder
# der Brunnen zu grob oder der Fahnenrand zu fein. Im Vorlauf (feines Band
# nur bis 70 m) war die Zelle bei r = 140 m - genau dort, wo die 1-K-Front
# steht - 14.7 m breit; die Front lag in einer einzigen Zelle.
C["mesh"].update({
    "n_r_well":     6,       # ueber den Brunnenradius (0.1 m)
    "r_near_m":    20.0,     # Nahfeld, stark gestaffelt: 0.06 -> 1.0 m
    "n_r_near":    60,
    "bias_r_near":  1.05,
    "r_fine_m":   180.0,     # Fahnenband, nahezu gleichmaessig 1.0 -> 2.0 m
    "n_r_fine":   110,       # (die 1-K-Front steht nach 30 a bei ~150 m)
    "bias_r_fine":  1.0065,
    "n_r_far":     40,       # Fernfeld 180 -> 1000 m
    "bias_r_far":   1.088,
    "n_z_aquifer": 32,                 # 38 m / 32 = 1.19 m
    "caprock_fine_depth_m": 80.0,      # feine Zone deckt die Leitfront ab
    "n_z_caprock_fine":     40,        # 2.0 m
    "n_z_caprock_coarse":   12,
})

# --- Zeit und Ausgabe -------------------------------------------------
C["time"]["dt_seconds"] = 86400.0
C["time"]["output_every_n_steps"] = 10
# Absolut neben diese Datei, nicht relativ zum Arbeitsverzeichnis:
# sonst landen die Ergebnisse dort, wo man das Skript zufaellig
# aufgerufen hat.
C["output"]["out_dir"] = str(Path(__file__).resolve().parent / "ergebnisse_2d")
C["output"]["prefix"] = "ates2d"


if __name__ == "__main__":
    ein = sum(p for p in P if p > 0)
    aus = sum(-p for p in P if p < 0)
    kf = _a["kf_m_s"]
    print("FALL 1 - ATES 2D radialsymmetrisch, OHNE Grundwasserstroemung")
    print(f"  Aquifer {_a['maechtigkeit_m']:.0f} m, kf = {kf:.2e} m/s "
          f"(k = {_k(kf):.3e} m2), Deckgestein 2 x {_d['maechtigkeit_m']:.0f} m")
    print(f"  Einspeisung {ein/1e3:7.1f} kW-Monate | "
          f"Entnahme {aus/1e3:7.1f} kW-Monate")
    print(f"  -> Deckungsgrad ist bei {100*ein/aus:.1f} % gedeckelt")
    print(f"  Referenz-Massenstrom "
          f"{C['operation']['mass_flow_rate_kg_s']:.2f} kg/s")
    print()
    sys.exit(main())

