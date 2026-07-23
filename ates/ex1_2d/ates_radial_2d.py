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
from __future__ import annotations

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
    h_screen = (cfg["layers"]["aquifer_thickness_m"]
                - cfg["well"]["screen_top_offset_m"]
                - cfg["well"]["screen_bottom_offset_m"])
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

    figdir = out_dir / "figures"; figdir.mkdir(exist_ok=True)
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
    z0 = max(0, times[-1] - 3*YEAR)
    plot_T(axs[1], z0, times[-1], False)
    axs[1].set_title(f"Zoom: letzte 3 Betriebsjahre")
    fig.tight_layout()
    fig.savefig(figdir / "T_vs_time.png", dpi=130); plt.close(fig)
    print("  saved T_vs_time.png")

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
            ax.set_xlabel("Betriebsjahr"); ax.set_ylabel("Recovery-Effizienz R [%]")
            ax.set_title("ATES 2D radial — jährliche Recovery-Effizienz\n"
                         f"R = (T̄_Förder − T_amb)/(T_inj − T_amb),  T_inj={T_inj-273.15:.0f} °C, "
                         f"T_amb={T_amb-273.15:.0f} °C")
            ax.set_ylim(0, max(100, max(R_years)*100*1.15)); ax.grid(alpha=0.3, axis="y")
            ax.axhline(100, color="k", lw=0.5, ls=":")
            fig.tight_layout()
            fig.savefig(figdir / "recovery_efficiency.png", dpi=130); plt.close(fig)
            print("  saved recovery_efficiency.png")
            print("  Recovery R pro Jahr [%]:", [f"{r*100:.0f}" for r in R_years])

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
            z0 = max(0, times[-1] - 4*YEAR)
            sel = times >= z0
            ax.plot(times[sel], demand_p[sel]/1e3, color="0.4", lw=1.6, ls="--",
                    label="Geforderte Förderleistung |P|")
            ax.plot(times[sel], deliver[sel]/1e3, color="#1f77b4", lw=1.6,
                    label="Gelieferte Förderleistung (bedarfsgeführt)")
            ax.set_xlabel("Zeit [Tage]"); ax.set_ylabel("Förderleistung [kW]")
            ax.set_title("ATES 2D radial — Bedarfsgeführte Förderung: gelieferte vs. geforderte Leistung "
                         "(letzte 4 Jahre)")
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
    return 0


if __name__ == "__main__":
    sys.exit(main())
