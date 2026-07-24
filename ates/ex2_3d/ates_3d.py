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
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

import gmsh
import re
import numpy as np

# ======================================================================
#  CONFIG  --  hier alles anpassen
# ======================================================================
CONFIG: dict = {
    # An MOOSE-ATES + das validierte 2D-Modell angelehnt, aber echtes 3D mit
    # regionaler Grundwasserströmung (die Fahne driftet stromab).
    "domain": {
        "size_x_m":   400.0,    # Strömungsrichtung (Platz für Fahnendrift)
        "size_y_m":   250.0,    # quer
        "z_base_m":     0.0,    # untere Modellgrenze (z-Koordinate)
    },
    "layers": {
        # Deckgestein dick genug, dass die Wärmeleitfront den auf T_amb fixierten
        # Rand nicht erreicht (bei 2 Jahren ~18 m -> 40 m je Seite reicht klar).
        "caprock_bottom_thickness_m": 40.0,
        "aquifer_thickness_m":        20.0,
        "caprock_top_thickness_m":    40.0,
    },
    "wells": {
        # Single-Well-Anlage. Brunnen als kleine Filtersäule über die volle
        # Aquiferhöhe (Offsets 0). Der Lateralrand des Aquifers ist Druck-Outlet
        # (bzw. GW-Gradient), damit injiziertes Wasser entweichen kann.
        "hot_well_xy":   ( 0.0,  0.0),     # (x, y) Lage des Brunnens
        "screen_bottom_offset_m": 0.0,     # Filter über volle Aquiferhöhe
        "screen_top_offset_m":    0.0,
        "screen_dx_m":             1.0,    # x-Ausdehnung des Filtervolumens
        "screen_dy_m":             1.0,    # y-Ausdehnung des Filtervolumens
        "screen_permeability_m2":  1.0e-11,# = Aquifer (Filter voll durchlässig)
        # Förderregelung wie im 2D: "fixed" (feste Rate) oder "demand"
        # (bedarfsgeführt, iterativ). Im 3D ist ein Lauf teuer (~1.5 h) →
        # "fixed" als Default; "demand" macht n Wiederholungsläufe.
        "production_control":      "fixed",
        "max_rate_factor":         6.0,
        "demand_iterations":       3,
    },
    # ------------------------------------------------------------------
    # REGIONALE GRUNDWASSERSTRÖMUNG  (das 3D-Kernmerkmal)
    # ------------------------------------------------------------------
    # Ein linearer Druckgradient auf der Lateral-Aquifer-Fläche prägt eine
    # großräumige Hintergrundströmung auf → die thermische Fahne driftet in
    # Strömungsrichtung. Kombiniert mit dem hydrostatischen Vertikaldruck.
    # ------------------------------------------------------------------
    "regional_gw": {
        "enable":            True,
        "gradient_m_per_m":  2.0e-3,       # hydraulischer Gradient (dimensionslos)
        "direction_deg":     0.0,          # 0° = +x, 90° = +y, ...
    },
    # Strukturiertes/feineres Netz: fein am Brunnen und stromab (Driftbahn),
    # gröber im Fernfeld.  (Tetraeder mit Distance-Verfeinerung, konforme
    # Schichtgrenzen durch OCC-Fragmentierung.)
    "mesh": {
        "size_in_well_m":       0.5,
        "size_near_wells_m":    1.7,       # Fahnen-/Driftregion fein
        "size_far_m":          20.0,       # Fernfeld grob (plumefern)
        "well_size_radius_m":  22.0,       # fein bis hierher (deckt Fahne + Drift)
        "well_size_radius_far_m": 65.0,    # ab dem: grob
    },
    "materials": {
        # MOOSE-Werte: Aquifer & Deckgestein gleiche thermische Eigenschaften,
        # nur Permeabilität unterscheidet sich (Deckgestein ~strömungsdicht).
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
        # Auftrieb wie im 2D: T-abhängige Dichte & Viskosität.
        "rho_ref_kg_m3":  1000.0,
        "T_ref_K":         283.15,
        "beta_1_per_K":   -4.0e-4,     # thermische Ausdehnung (heißes Wasser leichter)
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
        "T_K":  283.15,      # 10 °C Umgebung (= fluid.T_ref)
        "p_Pa": 1.0e5,       # Referenzdruck Oberkante (hydrostatisch nach unten)
    },
    "operation": {
        # T_hot_K = Injektionstemperatur T_inj; Massenstrom folgt aus der
        # Monatsleistung ṁ = P/(cp·(T_inj−T_amb)).
        "mass_flow_rate_kg_s": 3.0,    # nur Referenz im 4-Phasen-Modus
        "T_hot_K":  333.15,            # T_inj = 60 °C
        "T_cold_K": 283.15,
        "fluid_storage_1_per_Pa": 4.5e-10,
        "solid_storage_1_per_Pa": 1.0e-10,
    },
    # ------------------------------------------------------------------
    # ZYKLEN – HIER FÜR STUDIERENDE
    # ------------------------------------------------------------------
    # Ein vollständiger Zyklus besteht aus 4 aufeinander folgenden Phasen:
    #   1) charge                    – Beladung (heiß injizieren)
    #   2) storage_after_charge      – Pause/Speicherung nach Beladung
    #   3) discharge                 – Förderung (Wasser entnehmen)
    #   4) storage_after_discharge   – Pause/Speicherung nach Förderung
    #
    # Periode T_Zyklus = charge + storage_after_charge + discharge + storage_after_discharge
    # Gesamt­simulations­zeit = n_cycles * T_Zyklus
    #
    # Beispiele:
    #   - Saisonal (1 Jahr/Zyklus): 91.25 / 91.25 / 91.25 / 91.25
    #   - Sommerladung 120 d, Winterförderung 120 d, sonst Pause: 120 / 60 / 120 / 60
    #   - Phase auf 0 setzen, um sie zu deaktivieren.
    #
    # Alternativ — Modus B (Monatsprofil): Setze cycles.monthly_power_W auf eine
    # Liste von 12 Monatsleistungen [W] (positiv = laden/injizieren, negativ =
    # fördern/entnehmen, 0 = Stillstand). Dann wird die 4-Phasen-Logik
    # überschrieben; jeder Monat dauert 365.25/12 ≈ 30.44 d und die Sequenz wird
    # n_cycles-mal (= Jahre) wiederholt. Als Referenzleistung dient
    #   P_ref = mass_flow_rate_kg_s · cp · (T_hot_K − T_K),
    # d. h. der Massenstrom wird linear mit P_month/P_ref skaliert. Optional gibt
    # cycles.monthly_T_inj_K (12 Werte) die monatliche Vorlauftemperatur der
    # Beladung vor (sonst T_hot_K). Auf None lassen für Modus A.
    # ------------------------------------------------------------------
    "cycles": {
        "n_cycles":                     2,     # Betriebsjahre (Modus B)
        "charge_days":                 90,     # Modus A: Beladung (Tage)
        "storage_after_charge_days":    0,
        "discharge_days":              90,
        "storage_after_discharge_days": 0,
        "ramp_days":                    3.0,   # Übergangsrampe zwischen Monaten
        # --- Modus B: Monatsprofil (AKTIV) — P[W] Jan…Dez, ΣP≈0 ---
        # Aus P folgt die Pumprate ṁ = P/(cp·(T_inj−T_amb)); Spitze ~0.5 MW
        # → ṁ ≈ 2.4 kg/s. Beliebig editierbar.
        "monthly_power_W": [
            -400_000, -350_000, -120_000, 0,
            +250_000, +450_000, +500_000, +400_000, +120_000,
            -180_000, -320_000, -350_000,
        ],
        "monthly_T_inj_K":              None,   # optional: 12 Vorlauftemperaturen [K]
    },
    "time": {
        "dt_seconds":           86400.0,    # 1 Tag (feine Zeitauflösung)
        "output_every_n_steps": 5,
        "gravity":              True,        # Auftrieb (heißes Wasser steigt)
    },
    "output": {
        "prefix":    "ates_3d",
        "out_dir":   "out",
        "variables": ["T", "p", "darcy_velocity"],
    },
    "solver": {
        # BiCGSTAB+ILUT (schnell) oder "SparseLU" (direkt, robust).
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

# ======================================================================
#  Mesh – gmsh
# ======================================================================
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

    # Alles fragmentieren (konforme Schnittflächen)
    gmsh.model.occ.fragment([(3, box_cb)],
                            [(3, box_aq), (3, box_ct), (3, box_hw)])
    gmsh.model.occ.synchronize()

    # Volumen klassifizieren
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
            # innerhalb Aquifer-z – evtl. Brunnenbox?
            xc = 0.5 * (xmin + xmax)
            yc = 0.5 * (ymin + ymax)
            small = ext_x < 0.5 * Lx
            if small and abs(xc - hw[0]) < dx and abs(yc - hw[1]) < dy:
                vol_hw.append(tag)
            else:
                vol_aq.append(tag)
    if not vol_aq or not vol_hw:
        raise RuntimeError("Volumenklassifizierung fehlgeschlagen (gmsh-Fragmentierung).")

    # Top-/Bottom-Außenflächen + Lateral-Aquifer-Flächen
    surf_top, surf_bot, surf_lat_aq = [], [], []
    for dim, tag in gmsh.model.getEntities(2):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(dim, tag)
        zc = 0.5 * (zmin + zmax)
        # Top/Bottom: große horizontale Außenflächen
        if (xmax - xmin) >= 0.9 * Lx and abs(zc - z_top) < 1e-6:
            surf_top.append(tag)
            continue
        if (xmax - xmin) >= 0.9 * Lx and abs(zc - z_base) < 1e-6:
            surf_bot.append(tag)
            continue
        # Lateral-Aquifer: vertikale Außenflächen, deren z-Bereich genau die Aquiferdicke umfasst
        on_x_edge = abs(xmin - x0) < 1e-6 and abs(xmax - x0) < 1e-6
        on_x_edge_pos = abs(xmin - (x0 + Lx)) < 1e-6 and abs(xmax - (x0 + Lx)) < 1e-6
        on_y_edge = abs(ymin - y0) < 1e-6 and abs(ymax - y0) < 1e-6
        on_y_edge_pos = abs(ymin - (y0 + Ly)) < 1e-6 and abs(ymax - (y0 + Ly)) < 1e-6
        on_outer = on_x_edge or on_x_edge_pos or on_y_edge or on_y_edge_pos
        in_aquifer_z = (zmin >= z_aq_bot - 1e-6) and (zmax <= z_aq_top + 1e-6)
        if on_outer and in_aquifer_z:
            surf_lat_aq.append(tag)

    # Hüllflächen der Brunnenbox (für Distanzfeld / Neumann‑Wärme­fluss)
    surf_hw = []
    for tag in vol_hw:
        for d, t in gmsh.model.getBoundary([(3, tag)], oriented=False):
            if d == 2:
                surf_hw.append(abs(t))
    surf_hw = sorted(set(surf_hw))

    # Physical Groups (Reihenfolge -> tag -> MaterialID nach reindex)
    gmsh.model.addPhysicalGroup(3, vol_aq, tag=1, name="aquifer")
    gmsh.model.addPhysicalGroup(3, vol_ct, tag=2, name="caprock_top")
    gmsh.model.addPhysicalGroup(3, vol_cb, tag=3, name="caprock_bottom")
    gmsh.model.addPhysicalGroup(3, vol_hw, tag=4, name="hot_well_vol")
    gmsh.model.addPhysicalGroup(2, surf_top, tag=10, name="top")
    gmsh.model.addPhysicalGroup(2, surf_bot, tag=11, name="bottom")
    gmsh.model.addPhysicalGroup(2, surf_hw,  tag=12, name="hot_well_surf")
    if surf_lat_aq:
        gmsh.model.addPhysicalGroup(2, surf_lat_aq, tag=14, name="lateral_aquifer")

    # Hülle der Brunnenbox (für Distanzfeld)
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

    # Innerhalb der Brunnenbox sehr feines Netz (Punkte der Filterbox)
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
    return {
        "domain":          f"{prefix}_domain.vtu",
        "top":             f"{prefix}_physical_group_top.vtu",
        "bottom":          f"{prefix}_physical_group_bottom.vtu",
        "hot_well_vol":    f"{prefix}_physical_group_hot_well_vol.vtu",
        "hot_well_surf":   f"{prefix}_physical_group_hot_well_surf.vtu",
        "lateral_aquifer": f"{prefix}_physical_group_lateral_aquifer.vtu",
    }


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


# ======================================================================
#  Zyklus-Kurven
# ======================================================================
def build_cycle_curves(cfg: dict, rate_mult=None) -> dict:
    """Normierte Pump-/Leistungskurve g(t) ∈ [−1,1] + Beladungs-Intervalle
    (wie im validierten 2D-Modell).

    Brunnenmodell (physikalisch korrekt für OGS-HT): Massenquelle IMMER aktiv
    (∝ P, /ρ in build_prj); Injektionstemperatur per Dirichlet nur in den
    Beladungs-Intervallen; Förderung = reine Massensenke (dynamische Entnahme-
    temperatur). rate_mult skaliert die Förderrate je Monat (Bedarfsführung).
    """
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

    # === 4-Phasen-Modus ===
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


# ======================================================================
#  XML / .prj Generierung
# ======================================================================
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

    # Dichte: T-abhängig ρ(T)=ρ_ref(1+β(T−T_ref)) für Auftrieb (β<0), sonst const
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

    # Viskosität: optional T-abhängig (Linear) — verstärkt den Auftrieb
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

    # Effektive Wärmeleitfähigkeit: Mischung aus Phasen über Porosität
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

    # Brunnenfilter-Geometrie
    h_screen = (cfg["layers"]["aquifer_thickness_m"]
                - cfg["wells"]["screen_top_offset_m"]
                - cfg["wells"]["screen_bottom_offset_m"])
    dx_w = cfg["wells"]["screen_dx_m"]
    dy_w = cfg["wells"]["screen_dy_m"]
    V_well = dx_w * dy_w * h_screen

    # Basis-Amplitude der Massenquelle. WICHTIG: Der HT-Druck-Quellterm ist
    # VOLUMETRISCH [1/s], nicht massenbasiert -> durch ρ_f teilen (sonst ~1000×
    # Über-Injektion). cycle_power ∈[−1,1] skaliert zeitabhängig.
    rho_f = fluid["rho_ref_kg_m3"]
    q_mass_amp = curves["mdot_nom_kg_s"] / rho_f / V_well   # 1/s (= m³/m³/s)

    root = ET.Element("OpenGeoSysProject")

    # -- Meshes
    meshes = ET.SubElement(root, "meshes")
    _mesh_keys = ["domain", "top", "bottom", "hot_well_vol", "hot_well_surf"]
    if "lateral_aquifer" in mesh_files:
        _mesh_keys.append("lateral_aquifer")
    for k in _mesh_keys:
        _se(meshes, "mesh", mesh_files[k])

    # -- Processes
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

    # -- Media (Reihenfolge entspricht den MaterialIDs nach msh2vtu reindex)
    well_mat = dict(cfg["materials"]["aquifer"])
    well_mat["permeability_m2"] = cfg["wells"]["screen_permeability_m2"]
    disp = cfg["dispersion"]
    media = _se(root, "media")
    _add_medium(media, 0, cfg["materials"]["aquifer"],        fluid, op, disp)
    _add_medium(media, 1, cfg["materials"]["caprock_top"],    fluid, op, disp)
    _add_medium(media, 2, cfg["materials"]["caprock_bottom"], fluid, op, disp)
    _add_medium(media, 3, well_mat,                           fluid, op, disp)

    # -- Time loop
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

    # Adaptive Zeitschrittweite: IterationNumberBasedTimeStepping kann die
    # Schrittweite bei Nichtkonvergenz reduzieren (und danach wieder
    # erhöhen) — anders als FixedTimeStepping, das an harten Phasen-
    # übergängen (Flussumkehr + Dirichlet-T-Sprung am Brunnen) sonst
    # abbricht ("Time stepper cannot reduce the time step size further").
    # Normalbetrieb bleibt beim 1-Tages-Schritt (maximum_dt = dt0).
    dt0 = cfg["time"]["dt_seconds"]
    ts = _se(p_ref, "time_stepping")
    _se(ts, "type",        "IterationNumberBasedTimeStepping")
    _se(ts, "t_initial",   0.0)
    _se(ts, "t_end",       curves["t_total"])
    _se(ts, "initial_dt",  dt0)
    _se(ts, "minimum_dt",  dt0 / 64.0)   # erlaubt Reduktion an den Übergängen
    _se(ts, "maximum_dt",  dt0)          # nie gröber als der Standard-Schritt
    # Bei wenigen Nichtlinear-Iterationen wächst dt (bis maximum_dt), bei
    # vielen schrumpft es; bei Nichtkonvergenz wird der Schritt verworfen
    # und mit kleinerem dt wiederholt.
    _se(ts, "number_iterations", "1 4 8 12")
    _se(ts, "multiplier",        "1.5 1.0 0.5 0.25")

    # -- Output: feste Ausgabezeitpunkte -> gleichmäßige Snapshots,
    #    unabhängig von der jetzt variablen Schrittweite.
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

    # -- Parameters
    params = _se(root, "parameters")
    _const_param(params, "T0", init["T_K"])

    # Druck: bei Gravitation hydrostatisch p(z)=p_top+ρ_ref·g·(z_top−z); die
    # regionale GW-Strömung ist ein zusätzlicher lateraler Gradient in x/y.
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

    # Massenquelle (volumetrisch, /ρ) — Kurve cycle_power ∈[−1,1]
    _const_param(params, "q_mass_amp", q_mass_amp)
    _curve_scaled_param(params, "q_mass_well", "cycle_power", "q_mass_amp")
    # Injektionstemperatur(en) für die Beladungs-Dirichlet-BCs
    tinj_values = sorted({round(iv[2], 6) for iv in curves["charge_intervals"]})
    tinj_param = {}
    for idx, Ti in enumerate(tinj_values):
        nm = "T_inj" if len(tinj_values) == 1 else f"T_inj_{idx}"
        _const_param(params, nm, Ti); tinj_param[round(Ti, 6)] = nm

    # -- Curves
    cv = _se(root, "curves")
    t, v = curves["cycle_power"]; _curve_xml(cv, "cycle_power", t, v)

    # -- Process variables
    pvars = _se(root, "process_variables")

    # Temperatur
    pv_T = _se(pvars, "process_variable")
    _se(pv_T, "name", "T")
    _se(pv_T, "components", 1)
    _se(pv_T, "order", 1)
    _se(pv_T, "initial_condition", "T0")
    bcs_T = _se(pv_T, "boundary_conditions")
    # Aussenränder (Top/Bottom): konstante Hintergrund-T
    for face in ("top", "bottom"):
        bc = _se(bcs_T, "boundary_condition")
        _se(bc, "mesh",      Path(mesh_files[face]).stem)
        _se(bc, "type",      "Dirichlet")
        _se(bc, "parameter", "T0")
    # Brunnen: Injektionstemperatur T_inj NUR während der Beladungs-Intervalle
    # (DirichletWithinTimeInterval). Bei Förderung/Ruhe kein T-Zwang -> die
    # Entnahmetemperatur wird dynamisch berechnet.
    well_mesh_T = Path(mesh_files["hot_well_vol"]).stem
    for t0, t1, Ti in curves["charge_intervals"]:
        bc = _se(bcs_T, "boundary_condition")
        _se(bc, "mesh", well_mesh_T)
        _se(bc, "type", "DirichletWithinTimeInterval")
        _se(bc, "parameter", tinj_param[round(Ti, 6)])
        ti = _se(bc, "time_interval")
        _se(ti, "start", f"{t0:.6e}"); _se(ti, "end", f"{t1:.6e}")

    # Druck
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
    # Lateral-Aquifer als Druck-Outlet, damit das am Brunnen injizierte
    # Wasser entweichen kann. Falls regional_gw.enable=True: linearer
    # Druckgradient statt konstantem p0 → Hintergrund-Strömung.
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

    # -- Solvers
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
    # Solver konfigurierbar: iterativ (BiCGSTAB+ILUT, skaliert gut auf große Netze)
    # oder direkt (SparseLU, robust). scaling=true gleicht T/p-Größenordnungen aus.
    stype = sol.get("solver_type", "BiCGSTAB")
    _se(eig, "solver_type", stype)
    if stype != "SparseLU":
        _se(eig, "precon_type",        sol.get("precon_type", "ILUT"))
        _se(eig, "max_iteration_step", sol["linear_iter"])
        _se(eig, "error_tolerance",    sol["linear_tol"])
    _se(eig, "scaling",            "true")

    # Schreiben (mit XML-Deklaration und Einrückung)
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


# ======================================================================
#  OGS ausführen
# ======================================================================
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
    """|P|-gewichtete mittlere Brunnentemperatur je Fördermonat (P<0),
    gemittelt über die eingeschwungenen Jahre (zweite Laufhälfte)."""
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


# ======================================================================
#  CLI
# ======================================================================
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

    # Förderregelung: fest oder bedarfsgeführt (iterativ, wie 2D)
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
    return 0


if __name__ == "__main__":
    sys.exit(main())
