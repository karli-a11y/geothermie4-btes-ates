#!/usr/bin/env python3
"""
BTES 3D demo for OpenGeoSys 6 — Borehole Thermal Energy Storage.

Sondenfeld N x N im Untergrund. Je Sonde eine kleine Box als Subdomäne,
auf der eine volumetrische Wärmequelle aufgeprägt wird (positiv bei Beladung,
negativ bei Entladung). Wärmetransfer dominant über Wärmeleitung im Boden.

Im Gegensatz zum ATES (ates_3d.py) gibt es keinen Massenstrom — Strömungs-
permeabilität wird sehr niedrig gewählt, sodass die HT-Prozessgleichung
faktisch zu reiner Wärmeleitung degeneriert.

Aufruf:
    python btes_3d.py            # Mesh + .prj + OGS-Lauf
    python btes_3d.py --no-run   # nur Setup, kein OGS
    python btes_3d.py --no-mesh  # nur .prj erzeugen
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

import gmsh
import numpy as np

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
            fname = f"{output_prefix}_physical_group_{name}.vtu"
        mesh.save(str(output_path / fname), binary=True)


# ======================================================================
#  CONFIG  --  hier alles anpassen
# ======================================================================
CONFIG: dict = {
    "domain": {
        "size_x_m":  80.0,
        "size_y_m":  80.0,
        "z_base_m":   0.0,
    },
    "layers": {
        # Oberer "Cover" über dem Sondenfeld (z. B. Erdoberfläche)
        "soil_cover_thickness_m":   5.0,
        # Bereich, in dem die Sonden stecken (zwischen Cover und Bottom)
        "borehole_zone_thickness_m": 80.0,
        # Bottom unterhalb der Sonden (thermischer Puffer nach unten)
        "soil_bottom_thickness_m":  30.0,
    },
    "field": {
        # Sondenfeld – entweder N x N Raster (n_x, n_y, spacing_m) ODER
        # eine explizite Positionsliste positions=[(x,y), ...]
        "n_x":           3,
        "n_y":           3,
        "spacing_m":     5.0,
        "positions":     None,        # bei nicht-None überschreibt es das Raster
        # Tiefe je Sonde: top_offset bis bottom_offset (gemessen vom Top der borehole_zone)
        "top_offset_m":     2.0,
        "bottom_offset_m":  2.0,
        # Sonden als kleine Boxen modelliert (anstatt 1D Linien)
        "borehole_dx_m":    0.6,
        "borehole_dy_m":    0.6,
    },
    "mesh": {
        "size_in_borehole_m":    0.4,
        "size_near_field_m":     1.5,
        "size_far_m":           15.0,
        "field_size_radius_m":   8.0,    # bis hier feinmaschig
        "field_size_radius_far_m": 30.0,
    },
    "materials": {
        # Bodenmaterial (z. B. Lockergestein / Festgestein)
        "soil": {
            "permeability_m2":  1.0e-18,   # quasi impermeabel
            "porosity":          0.20,
            "rho_s_kg_m3":    2700.0,
            "cp_s_J_kgK":      900.0,
            "lambda_s_W_mK":     2.5,
        },
    },
    "fluid": {
        # Porenfluid (Wasser) — auch wenn k klein ist, brauchen wir die Phase
        "rho_ref_kg_m3":  1000.0,
        "T_ref_K":         293.15,
        "beta_1_per_K":    0.0,
        "viscosity_Pa_s":  1.0e-3,
        "cp_J_kgK":        4180.0,
        "lambda_W_mK":     0.6,
    },
    "dispersion": {
        "alpha_L_m": 0.0,    # ohne Strömung irrelevant
        "alpha_T_m": 0.0,
    },
    "initial": {
        "T_K":  283.15,
        "p_Pa": 0.0,
    },
    "operation": {
        # Heat-Source je Sonde: Leistung [W] (positiv = Beladung)
        "power_per_borehole_W":  2000.0,
        # Effektive Speicherzahlen (im HT-Prozess erforderlich)
        "fluid_storage_1_per_Pa": 4.5e-10,
        "solid_storage_1_per_Pa": 1.0e-10,
    },
    # ------------------------------------------------------------------
    # ZYKLEN – HIER FÜR STUDIERENDE
    # ------------------------------------------------------------------
    # Ein vollständiger Zyklus besteht aus 4 aufeinander folgenden Phasen:
    #   1) charge                    – Beladung (Wärme einspeisen)
    #   2) storage_after_charge      – Pause/Speicherung nach Beladung
    #   3) discharge                 – Förderung (Wärme entziehen)
    #   4) storage_after_discharge   – Pause/Speicherung nach Förderung
    #
    # Periode T_Zyklus = charge + storage_after_charge + discharge + storage_after_discharge
    # Gesamt­simulations­zeit = n_cycles * T_Zyklus
    #
    # Beispiele:
    #   - Saisonal (1 Jahr/Zyklus): 91.25 / 91.25 / 91.25 / 91.25
    #   - Sommerladung 120 d, Winterförderung 120 d, sonst Pause: 120 / 60 / 120 / 60
    #   - Phase auf 0 setzen, um sie zu deaktivieren.
    # ------------------------------------------------------------------
    "cycles": {
        "n_cycles":                        1,       # Anzahl Wiederholungen des Zyklus
        "charge_days":                     91.25,   # Phase 1: Beladung (Tage)
        "storage_after_charge_days":       91.25,   # Phase 2: Pause nach Beladung (Tage)
        "discharge_days":                  91.25,   # Phase 3: Förderung (Tage)
        "storage_after_discharge_days":    91.25,   # Phase 4: Pause nach Förderung (Tage)
        "ramp_days":                       7.0,     # Sanfte Übergangsrampe zwischen Phasen (Tage)
    },
    "time": {
        "dt_seconds":           7 * 86400.0,
        "output_every_n_steps": 1,
        "gravity":              False,
    },
    "output": {
        "prefix":    "btes_3d",
        "out_dir":   "out",
        "variables": ["T", "p", "darcy_velocity"],
    },
    "solver": {
        "linear_tol":      1.0e-12,
        "linear_iter":     10000,
        "nonlinear_iter":  20,
        "rel_tol_T":       1.0e-4,
        "rel_tol_p":       1.0e-4,
    },
}

DAY = 86400.0


# ======================================================================
#  Mesh — gmsh
# ======================================================================
def build_mesh(cfg: dict, out_dir: Path) -> Path:
    msh_path = out_dir / f"{cfg['output']['prefix']}.msh"

    Lx = cfg["domain"]["size_x_m"]
    Ly = cfg["domain"]["size_y_m"]
    z_base = cfg["domain"]["z_base_m"]

    t_bot   = cfg["layers"]["soil_bottom_thickness_m"]
    t_zone  = cfg["layers"]["borehole_zone_thickness_m"]
    t_cover = cfg["layers"]["soil_cover_thickness_m"]

    z_zone_bot = z_base + t_bot
    z_zone_top = z_zone_bot + t_zone
    z_top      = z_zone_top + t_cover

    dx_b = cfg["field"]["borehole_dx_m"]
    dy_b = cfg["field"]["borehole_dy_m"]
    z_bh_bot = z_zone_bot + cfg["field"]["bottom_offset_m"]
    z_bh_top = z_zone_top - cfg["field"]["top_offset_m"]
    h_bh = z_bh_top - z_bh_bot

    custom = cfg["field"].get("positions")
    if custom:
        bh_positions: list[tuple[float, float]] = [(float(x), float(y)) for x, y in custom]
    else:
        nx = cfg["field"]["n_x"]; ny = cfg["field"]["n_y"]; s = cfg["field"]["spacing_m"]
        bh_positions = []
        x0_f = -(nx - 1) * s / 2.0
        y0_f = -(ny - 1) * s / 2.0
        for ix in range(nx):
            for iy in range(ny):
                bh_positions.append((x0_f + ix * s, y0_f + iy * s))

    s_in   = cfg["mesh"]["size_in_borehole_m"]
    s_near = cfg["mesh"]["size_near_field_m"]
    s_far  = cfg["mesh"]["size_far_m"]
    r_near = cfg["mesh"]["field_size_radius_m"]
    r_far  = cfg["mesh"]["field_size_radius_far_m"]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("btes")

    x0, y0 = -Lx / 2.0, -Ly / 2.0
    box_bot   = gmsh.model.occ.addBox(x0, y0, z_base,     Lx, Ly, t_bot)
    box_zone  = gmsh.model.occ.addBox(x0, y0, z_zone_bot, Lx, Ly, t_zone)
    box_cover = gmsh.model.occ.addBox(x0, y0, z_zone_top, Lx, Ly, t_cover)

    bh_boxes = []
    for x, y in bh_positions:
        b = gmsh.model.occ.addBox(x - dx_b / 2.0, y - dy_b / 2.0, z_bh_bot,
                                  dx_b, dy_b, h_bh)
        bh_boxes.append(b)

    gmsh.model.occ.fragment(
        [(3, box_bot)],
        [(3, box_zone), (3, box_cover)] + [(3, b) for b in bh_boxes],
    )
    gmsh.model.occ.synchronize()

    vol_soil, vol_bh = [], []
    for dim, tag in gmsh.model.getEntities(3):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(dim, tag)
        ext = xmax - xmin
        zc = 0.5 * (zmin + zmax)
        small = ext < 0.3 * Lx
        if small and (z_bh_bot - 1e-3) <= zc <= (z_bh_top + 1e-3):
            vol_bh.append((tag, 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)))
        else:
            vol_soil.append(tag)
    if len(vol_bh) != len(bh_positions):
        raise RuntimeError(f"Erwartet {len(bh_positions)} Sonden, gefunden {len(vol_bh)}")

    # Sortiere Sonden-Volumen nach (x, y) — stabile Zuordnung zur Position
    vol_bh.sort(key=lambda t: (t[1], t[2]))
    vol_bh_tags = [t[0] for t in vol_bh]

    # Außenflächen
    surf_top, surf_bot, surf_lat_zone = [], [], []
    for dim, tag in gmsh.model.getEntities(2):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(dim, tag)
        zc = 0.5 * (zmin + zmax)
        if (xmax - xmin) >= 0.9 * Lx and abs(zc - z_top) < 1e-6:
            surf_top.append(tag)
            continue
        if (xmax - xmin) >= 0.9 * Lx and abs(zc - z_base) < 1e-6:
            surf_bot.append(tag)
            continue
        # Lateral der "borehole_zone"
        on_outer = (abs(xmin - x0) < 1e-6 and abs(xmax - x0) < 1e-6) \
                   or (abs(xmin - (x0 + Lx)) < 1e-6 and abs(xmax - (x0 + Lx)) < 1e-6) \
                   or (abs(ymin - y0) < 1e-6 and abs(ymax - y0) < 1e-6) \
                   or (abs(ymin - (y0 + Ly)) < 1e-6 and abs(ymax - (y0 + Ly)) < 1e-6)
        in_zone_z = (zmin >= z_zone_bot - 1e-6) and (zmax <= z_zone_top + 1e-6)
        if on_outer and in_zone_z:
            surf_lat_zone.append(tag)

    # Physical groups – fortlaufende Tags für Volumen, damit msh2vtu reindex
    # zusammenhängende MaterialIDs 0..N produziert
    gmsh.model.addPhysicalGroup(3, vol_soil, tag=1, name="soil")
    for i, tag in enumerate(vol_bh_tags):
        gmsh.model.addPhysicalGroup(3, [tag], tag=2 + i, name=f"bh_{i:02d}")
    gmsh.model.addPhysicalGroup(2, surf_top,       tag=100, name="top")
    gmsh.model.addPhysicalGroup(2, surf_bot,       tag=101, name="bottom")
    if surf_lat_zone:
        gmsh.model.addPhysicalGroup(2, surf_lat_zone, tag=102, name="lateral_zone")

    # Mesh-Größenfeld: feiner um Sondenfeld
    bh_surfaces = []
    for tag in vol_bh_tags:
        for d, t in gmsh.model.getBoundary([(3, tag)], oriented=False):
            if d == 2:
                bh_surfaces.append(abs(t))
    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "SurfacesList", list(set(bh_surfaces)))
    f_thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_thr, "InField", f_dist)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMin", s_near)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMax", s_far)
    gmsh.model.mesh.field.setNumber(f_thr, "DistMin", r_near)
    gmsh.model.mesh.field.setNumber(f_thr, "DistMax", r_far)
    gmsh.model.mesh.field.setAsBackgroundMesh(f_thr)

    bh_points = []
    for tag in vol_bh_tags:
        for d, t in gmsh.model.getBoundary([(3, tag)], recursive=True, oriented=False):
            if d == 0:
                bh_points.append((d, t))
    if bh_points:
        gmsh.model.mesh.setSize(bh_points, s_in)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

    gmsh.model.mesh.generate(3)
    gmsh.write(str(msh_path))
    gmsh.finalize()
    return msh_path


def _n_boreholes(cfg: dict) -> int:
    custom = cfg["field"].get("positions")
    return len(custom) if custom else cfg["field"]["n_x"] * cfg["field"]["n_y"]


def convert_mesh(cfg: dict, msh_path: Path, out_dir: Path) -> dict:
    pass  # msh2vtu provided at module level
    prefix = cfg["output"]["prefix"]
    msh2vtu(filename=msh_path, output_path=out_dir, output_prefix=prefix,
            dim=3, reindex=True, log_level="WARNING")
    n_bh = _n_boreholes(cfg)
    bh_meshes = {f"bh_{i:02d}": f"{prefix}_physical_group_bh_{i:02d}.vtu" for i in range(n_bh)}
    return {
        "domain":       f"{prefix}_domain.vtu",
        "top":          f"{prefix}_physical_group_top.vtu",
        "bottom":       f"{prefix}_physical_group_bottom.vtu",
        "lateral_zone": f"{prefix}_physical_group_lateral_zone.vtu",
        **bh_meshes,
        "_n_bh":        n_bh,
    }


# ======================================================================
#  Zyklen-Kurve
# ======================================================================
def build_cycle_curves(cfg: dict) -> dict:
    """Eine Kurve für alle Sonden: +1 Beladung, 0 Pause, -1 Förderung."""
    n     = cfg["cycles"]["n_cycles"]
    t_c   = cfg["cycles"]["charge_days"]                  * DAY
    t_sc  = cfg["cycles"]["storage_after_charge_days"]    * DAY
    t_d   = cfg["cycles"]["discharge_days"]               * DAY
    t_sd  = cfg["cycles"]["storage_after_discharge_days"] * DAY
    ramp  = max(60.0, cfg["cycles"]["ramp_days"] * DAY)

    # Phasen: (Dauer, q-Skalierung)
    phases = [
        ("charge",          t_c,  +1.0),
        ("storage_after_c", t_sc,  0.0),
        ("discharge",       t_d,  -1.0),
        ("storage_after_d", t_sd,  0.0),
    ]

    times = [0.0]
    vals  = [0.0]
    t_now = 0.0
    for _ in range(n):
        for _name, dur, q in phases:
            if dur <= 0.0: continue
            t_now += ramp
            times.append(t_now); vals.append(q)
            hold = max(0.0, dur - ramp)
            if hold > 0.0:
                t_now += hold
                times.append(t_now); vals.append(q)
    t_now += ramp
    times.append(t_now); vals.append(0.0)
    return {
        "t_total":  t_now,
        "cycle_q":  (np.array(times), np.array(vals)),
    }


# ======================================================================
#  XML / .prj
# ======================================================================
def _se(parent, tag, text=None, **attrs):
    el = ET.SubElement(parent, tag, **{k: str(v) for k, v in attrs.items()})
    if text is not None:
        el.text = str(text)
    return el


def _const_param(parent, name, value):
    p = _se(parent, "parameter")
    _se(p, "name", name)
    _se(p, "type", "Constant")
    _se(p, "value", value)


def _curve_scaled_param(parent, name, curve_name, base_param_name):
    p = _se(parent, "parameter")
    _se(p, "name", name)
    _se(p, "type", "CurveScaled")
    _se(p, "curve", curve_name)
    _se(p, "parameter", base_param_name)


def _add_const_property(parent, name, value):
    p = _se(parent, "property")
    _se(p, "name", name)
    _se(p, "type", "Constant")
    _se(p, "value", value)


def _add_phase_aqueous(phases, fluid, op):
    ph = _se(phases, "phase")
    _se(ph, "type", "AqueousLiquid")
    props = _se(ph, "properties")
    _add_const_property(props, "density",              fluid["rho_ref_kg_m3"])
    _add_const_property(props, "viscosity",            fluid["viscosity_Pa_s"])
    _add_const_property(props, "specific_heat_capacity", fluid["cp_J_kgK"])
    _add_const_property(props, "thermal_conductivity", fluid["lambda_W_mK"])
    _add_const_property(props, "storage",              op["fluid_storage_1_per_Pa"])


def _add_phase_solid(phases, mat, op):
    ph = _se(phases, "phase")
    _se(ph, "type", "Solid")
    props = _se(ph, "properties")
    _add_const_property(props, "density",                mat["rho_s_kg_m3"])
    _add_const_property(props, "specific_heat_capacity", mat["cp_s_J_kgK"])
    _add_const_property(props, "thermal_conductivity",   mat["lambda_s_W_mK"])
    _add_const_property(props, "storage",                op["solid_storage_1_per_Pa"])


def _add_medium(media, mid, mat, fluid, op, disp):
    med = _se(media, "medium", id=mid)
    phases = _se(med, "phases")
    _add_phase_aqueous(phases, fluid, op)
    _add_phase_solid(phases, mat, op)
    props = _se(med, "properties")
    _add_const_property(props, "porosity",     mat["porosity"])
    _add_const_property(props, "permeability", mat["permeability_m2"])
    p = _se(props, "property"); _se(p, "name", "thermal_conductivity")
    _se(p, "type", "EffectiveThermalConductivityPorosityMixing")
    _add_const_property(props, "thermal_longitudinal_dispersivity", disp["alpha_L_m"])
    _add_const_property(props, "thermal_transversal_dispersivity",  disp["alpha_T_m"])
    _add_const_property(props, "storage", 0.0)


def _curve_xml(parent, name, t, v):
    c = _se(parent, "curve")
    _se(c, "name", name)
    _se(c, "coords", " ".join(f"{x:.6e}" for x in t))
    _se(c, "values", " ".join(f"{x:.6e}" for x in v))


def _indent(elem, level=0):
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


def build_prj(cfg: dict, out_dir: Path, mesh_files: dict, curves: dict) -> Path:
    prefix = cfg["output"]["prefix"]
    fluid  = cfg["fluid"]
    op     = cfg["operation"]
    init   = cfg["initial"]
    sol    = cfg["solver"]

    # Wärmequelle pro Volumen [W/m³] je Sonde
    h_bh = (cfg["layers"]["borehole_zone_thickness_m"]
            - cfg["field"]["top_offset_m"] - cfg["field"]["bottom_offset_m"])
    V_bh = cfg["field"]["borehole_dx_m"] * cfg["field"]["borehole_dy_m"] * h_bh
    q_v = op["power_per_borehole_W"] / V_bh

    root = ET.Element("OpenGeoSysProject")

    meshes = _se(root, "meshes")
    _se(meshes, "mesh", mesh_files["domain"])
    _se(meshes, "mesh", mesh_files["top"])
    _se(meshes, "mesh", mesh_files["bottom"])
    if mesh_files.get("lateral_zone"):
        _se(meshes, "mesh", mesh_files["lateral_zone"])
    for i in range(mesh_files["_n_bh"]):
        _se(meshes, "mesh", mesh_files[f"bh_{i:02d}"])

    # Process
    processes = _se(root, "processes")
    proc = _se(processes, "process")
    _se(proc, "name", "HT"); _se(proc, "type", "HT"); _se(proc, "integration_order", 2)
    pv = _se(proc, "process_variables"); _se(pv, "temperature", "T"); _se(pv, "pressure", "p")
    sv = _se(proc, "secondary_variables")
    _se(sv, "secondary_variable", internal_name="darcy_velocity", output_name="darcy_velocity")
    _se(proc, "specific_body_force", "0 0 0")

    # Media: alle Subdomänen erhalten gleiches Soil-Material
    # MaterialIDs nach msh2vtu reindex: soil=0, bh_00=1, bh_01=2, ...
    media = _se(root, "media")
    soil = cfg["materials"]["soil"]
    disp = cfg["dispersion"]
    _add_medium(media, 0, soil, fluid, op, disp)
    for i in range(mesh_files["_n_bh"]):
        _add_medium(media, i + 1, soil, fluid, op, disp)

    # Time loop
    tl = _se(root, "time_loop")
    procs = _se(tl, "processes")
    p_ref = _se(procs, "process", ref="HT")
    _se(p_ref, "nonlinear_solver", "basic_picard")
    cc = _se(p_ref, "convergence_criterion")
    _se(cc, "type", "PerComponentDeltaX")
    _se(cc, "norm_type", "NORM2")
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

    # Parameters
    params = _se(root, "parameters")
    _const_param(params, "T0",   init["T_K"])
    _const_param(params, "p0",   init["p_Pa"])
    _const_param(params, "q_v_amp", q_v)
    _curve_scaled_param(params, "q_v_borehole", "cycle_q", "q_v_amp")

    cv = _se(root, "curves")
    t, v = curves["cycle_q"]
    _curve_xml(cv, "cycle_q", t, v)

    # Process variables
    pvars = _se(root, "process_variables")

    # T
    pv_T = _se(pvars, "process_variable")
    _se(pv_T, "name", "T"); _se(pv_T, "components", 1); _se(pv_T, "order", 1)
    _se(pv_T, "initial_condition", "T0")
    bcs = _se(pv_T, "boundary_conditions")
    for face in ("top", "bottom"):
        bc = _se(bcs, "boundary_condition")
        _se(bc, "mesh", Path(mesh_files[face]).stem)
        _se(bc, "type", "Dirichlet")
        _se(bc, "parameter", "T0")
    sts = _se(pv_T, "source_terms")
    for i in range(mesh_files["_n_bh"]):
        st = _se(sts, "source_term")
        _se(st, "mesh", Path(mesh_files[f"bh_{i:02d}"]).stem)
        _se(st, "type", "Volumetric")
        _se(st, "parameter", "q_v_borehole")

    # p (statisch, nur Dirichlet Boundary)
    pv_p = _se(pvars, "process_variable")
    _se(pv_p, "name", "p"); _se(pv_p, "components", 1); _se(pv_p, "order", 1)
    _se(pv_p, "initial_condition", "p0")
    bcs = _se(pv_p, "boundary_conditions")
    for face in ("top", "bottom"):
        bc = _se(bcs, "boundary_condition")
        _se(bc, "mesh", Path(mesh_files[face]).stem)
        _se(bc, "type", "Dirichlet")
        _se(bc, "parameter", "p0")

    # Solvers
    nls = _se(root, "nonlinear_solvers")
    n = _se(nls, "nonlinear_solver")
    _se(n, "name", "basic_picard"); _se(n, "type", "Picard")
    _se(n, "max_iter", sol["nonlinear_iter"]); _se(n, "linear_solver", "general_linear_solver")

    lss = _se(root, "linear_solvers")
    ls = _se(lss, "linear_solver")
    _se(ls, "name", "general_linear_solver")
    eig = _se(ls, "eigen")
    _se(eig, "solver_type", "BiCGSTAB"); _se(eig, "precon_type", "ILUT")
    _se(eig, "max_iteration_step", sol["linear_iter"])
    _se(eig, "error_tolerance",    sol["linear_tol"])
    _se(eig, "scaling", "true")

    _indent(root)
    prj_path = out_dir / f"{prefix}.prj"
    ET.ElementTree(root).write(prj_path, encoding="ISO-8859-1", xml_declaration=True)
    return prj_path


# ======================================================================
#  Run + CLI
# ======================================================================
def run_ogs(prj_path: Path) -> int:
    ogs_exe = shutil.which("ogs") or shutil.which("ogs.exe")
    if not ogs_exe:
        print("ogs.exe nicht im PATH", file=sys.stderr); return 1
    cmd = [ogs_exe, str(prj_path), "-o", str(prj_path.parent)]
    print(">>", " ".join(cmd))
    return subprocess.call(cmd)


def main() -> int:
    ap = argparse.ArgumentParser(description="BTES 3D OGS demo")
    ap.add_argument("--no-mesh", action="store_true")
    ap.add_argument("--no-run",  action="store_true")
    args = ap.parse_args()

    out_dir = Path(CONFIG["output"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = CONFIG["output"]["prefix"]
    msh_path = out_dir / f"{prefix}.msh"

    if not args.no_mesh:
        print("[1/3] gmsh: Sondenfeld ...")
        build_mesh(CONFIG, out_dir)
        print(f"      {msh_path}")
        print("[2/3] msh2vtu: Konvertierung ...")
        mesh_files = convert_mesh(CONFIG, msh_path, out_dir)
    else:
        n_bh = _n_boreholes(CONFIG)
        mesh_files = {
            "domain":       f"{prefix}_domain.vtu",
            "top":          f"{prefix}_physical_group_top.vtu",
            "bottom":       f"{prefix}_physical_group_bottom.vtu",
            "lateral_zone": f"{prefix}_physical_group_lateral_zone.vtu",
            "_n_bh":        n_bh,
            **{f"bh_{i:02d}": f"{prefix}_physical_group_bh_{i:02d}.vtu" for i in range(n_bh)},
        }

    print("[3/3] OGS-Projektdatei ...")
    curves = build_cycle_curves(CONFIG)
    prj_path = build_prj(CONFIG, out_dir, mesh_files, curves)
    print(f"      {prj_path}  (t_end = {curves['t_total']/DAY:.1f} d, {mesh_files['_n_bh']} Sonden)")

    if args.no_run:
        return 0
    print(">>> OGS starten")
    return run_ogs(prj_path)


if __name__ == "__main__":
    sys.exit(main())
