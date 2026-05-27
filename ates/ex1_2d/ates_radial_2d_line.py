#!/usr/bin/env python3
# coding: utf-8
"""
ATES 2D radialsymmetrisch — Variante mit LINIEN-QUELLE.

STATUS: SKIZZE. Mesh- und PRJ-Generierung sind implementiert, ein
end-to-end-OGS-Lauf ist noch nicht validiert.

Unterschied zur Variante `ates_radial_2d.py` (Volumen­quelle):

- Die Filterstrecke wird als **1D-Linie** auf der Symmetrie-Achse
  *r* = 0 modelliert, nicht als 2D-Rechteck mit Volumen­quell­term.
- Massen­injektion erfolgt als `NodalSourceTerm` **am obersten Knoten**
  dieser Linie. Innerhalb der Filterstrecke verteilt sich das Wasser
  numerisch über die Druckgleichung und tritt seitlich aus.
- Die Temperatur am Top-Knoten wird per Dirichlet-T-BC auf
  *T*<sub>inj</sub> gehalten — dadurch entspricht das Modell näher
  einer realen Brunnen­installation, bei der Wasser oben einläuft
  und durch die Filterstrecke nach außen strömt.

OGS-seitig wird der Top-Knoten als eigene Physical-Group (Point)
exportiert; OGS akzeptiert dort `NodalSourceTerm` für Massen- und
Energie-Eintrag.
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


def msh2vtu(filename, output_path, output_prefix, dim, reindex=True, log_level="WARNING"):
    import ogstools as ot
    from pathlib import Path as _P
    meshes = ot.Meshes.from_gmsh(filename=str(filename), dim=dim,
                                  reindex=reindex, log=False)
    output_path = _P(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    for name, mesh in meshes.items():
        fname = (f"{output_prefix}_domain.vtu" if name == "domain"
                 else f"{output_prefix}_physical_group_{name}.vtu")
        mesh.save(str(output_path / fname), binary=True)


# ======================================================================
#  CONFIG
# ======================================================================
CONFIG: dict = {
    "domain": {"r_max_m": 500.0},
    "layers": {
        "caprock_bottom_thickness_m": 60.0,
        "aquifer_thickness_m":        30.0,
        "caprock_top_thickness_m":    60.0,
    },
    "well": {
        # Filterstrecke als 1D-Linie auf der Achse r=0
        # — Studis konfigurieren Tiefen­bereich (gemessen von OK)
        "screen_top_depth_m":     65.0,   # Tiefe Filter-Oberkante
        "screen_bottom_depth_m":  85.0,   # Tiefe Filter-Unterkante
        # Im Linien-Modell hat das Bohrloch nominell keinen Radius;
        # für die Verfeinerung um die Achse verwenden wir trotzdem
        # einen Verfeinerungs­radius.
        "mesh_refine_radius_m":   3.0,
    },
    "mesh": {
        "size_in_filter_m":       0.4,
        "size_near_well_m":       1.0,
        "size_far_m":            15.0,
        "well_size_radius_m":     5.0,
        "well_size_radius_far_m": 120.0,
    },
    "materials": {
        "aquifer":         {"permeability_m2": 1.0e-12, "porosity": 0.25,
                            "rho_s_kg_m3": 2650.0, "cp_s_J_kgK": 1000.0,
                            "lambda_s_W_mK": 3.0},
        "caprock_top":     {"permeability_m2": 1.0e-18, "porosity": 0.05,
                            "rho_s_kg_m3": 2700.0, "cp_s_J_kgK":  900.0,
                            "lambda_s_W_mK": 2.0},
        "caprock_bottom":  {"permeability_m2": 1.0e-18, "porosity": 0.05,
                            "rho_s_kg_m3": 2700.0, "cp_s_J_kgK":  900.0,
                            "lambda_s_W_mK": 2.0},
    },
    "fluid": {
        "rho_ref_kg_m3":  1000.0,
        "T_ref_K":         293.15,
        "beta_1_per_K":    0.0,
        "viscosity_Pa_s":  1.0e-3,
        "cp_J_kgK":        4180.0,
        "lambda_W_mK":     0.6,
    },
    "dispersion": {"alpha_L_m": 5.0, "alpha_T_m": 1.0},
    "initial":    {"T_K": 283.15, "p_Pa": 0.0,
                   "T_surface_K": 283.15,
                   "geothermal_gradient_K_per_m": 0.0},
    "operation": {
        # Massenstrom am Top-Knoten der Filterstrecke
        "mass_flow_rate_kg_s": 0.5,
        "T_hot_K":  353.15,
        "T_cold_K": 283.15,
        "fluid_storage_1_per_Pa": 4.5e-10,
        "solid_storage_1_per_Pa": 1.0e-10,
    },
    "cycles": {
        "n_cycles":                       1,
        "charge_days":                   91.25,
        "storage_after_charge_days":     91.25,
        "discharge_days":                91.25,
        "storage_after_discharge_days":  91.25,
        "ramp_days":                      7.0,
        "monthly_power_W":   None,
        "monthly_T_inj_K":   None,
    },
    "time":   {"dt_seconds": 7*86400.0, "output_every_n_steps": 1, "gravity": False},
    "output": {"prefix": "ates_radial_2d_line", "out_dir": "out",
               "variables": ["T", "p", "darcy_velocity"]},
    "solver": {"linear_tol": 1e-12, "linear_iter": 10000,
               "nonlinear_iter": 20, "rel_tol_T": 1e-4, "rel_tol_p": 1e-4},
}

DAY = 86400.0


# ======================================================================
#  Mesh — gmsh (2D radial mit eingebetteter 1D-Filter-Linie)
# ======================================================================
def build_mesh(cfg: dict, out_dir: Path) -> Path:
    msh_path = out_dir / f"{cfg['output']['prefix']}.msh"
    r_max = cfg["domain"]["r_max_m"]
    t_cb  = cfg["layers"]["caprock_bottom_thickness_m"]
    t_aq  = cfg["layers"]["aquifer_thickness_m"]
    t_ct  = cfg["layers"]["caprock_top_thickness_m"]
    z_top = t_cb + t_aq + t_ct

    s_top_d = cfg["well"]["screen_top_depth_m"]
    s_bot_d = cfg["well"]["screen_bottom_depth_m"]
    z_filter_top = z_top - s_top_d   # interne z-Koordinate
    z_filter_bot = z_top - s_bot_d
    if not (0 <= z_filter_bot < z_filter_top <= z_top):
        raise ValueError("Filterstrecke außerhalb der Schicht­geometrie.")
    if not (t_cb <= z_filter_bot and z_filter_top <= t_cb + t_aq):
        print("WARNUNG: Filterstrecke überschreitet die Aquifer-Schicht.")

    m = cfg["mesh"]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("ates_line")

    r_cb = gmsh.model.occ.addRectangle(0, 0,    0, r_max, t_cb)
    r_aq = gmsh.model.occ.addRectangle(0, t_cb, 0, r_max, t_aq)
    r_ct = gmsh.model.occ.addRectangle(0, t_cb+t_aq, 0, r_max, t_ct)

    # Filter als 1D-Linie auf r=0
    p_top = gmsh.model.occ.addPoint(0, z_filter_top, 0)
    p_bot = gmsh.model.occ.addPoint(0, z_filter_bot, 0)
    filter_line = gmsh.model.occ.addLine(p_top, p_bot)

    gmsh.model.occ.fragment([(2, r_cb)], [(2, r_aq), (2, r_ct), (1, filter_line)])
    gmsh.model.occ.synchronize()

    # Klassifikation der 2D-Flächen nach z-Mittelpunkt
    surf_cb, surf_aq, surf_ct = [], [], []
    for dim, tag in gmsh.model.getEntities(2):
        xmin, ymin, _, xmax, ymax, _ = gmsh.model.occ.getBoundingBox(dim, tag)
        zc = 0.5*(ymin + ymax)
        if   zc < t_cb:                 surf_cb.append(tag)
        elif zc < t_cb + t_aq:          surf_aq.append(tag)
        else:                           surf_ct.append(tag)

    # Filter-Linie (kann durch fragment in Stücke zerteilt sein, alle einsammeln)
    line_filter = []
    for dim, tag in gmsh.model.getEntities(1):
        xmin, ymin, _, xmax, ymax, _ = gmsh.model.occ.getBoundingBox(dim, tag)
        if abs(xmin) < 1e-6 and abs(xmax) < 1e-6:
            if (z_filter_bot - 1e-6) <= 0.5*(ymin+ymax) <= (z_filter_top + 1e-6):
                line_filter.append(tag)

    # Top-Punkt der Filter-Linie (höchster z-Wert auf r=0)
    top_point = None
    best_z = -1e9
    for dim, tag in gmsh.model.getEntities(0):
        xmin, ymin, _, *_ = gmsh.model.occ.getBoundingBox(dim, tag)
        if abs(xmin) < 1e-6 and abs(ymin - z_filter_top) < 1e-3:
            if ymin > best_z:
                best_z = ymin; top_point = tag
    if top_point is None:
        raise RuntimeError("Filter-Top-Punkt nicht gefunden.")

    # Randkanten
    edge_top, edge_bot, edge_far, edge_axis = [], [], [], []
    for dim, tag in gmsh.model.getEntities(1):
        xmin, ymin, _, xmax, ymax, _ = gmsh.model.occ.getBoundingBox(dim, tag)
        if abs(ymin - z_top) < 1e-6 and abs(ymax - z_top) < 1e-6:
            edge_top.append(tag)
        elif abs(ymin) < 1e-6 and abs(ymax) < 1e-6:
            edge_bot.append(tag)
        elif abs(xmin - r_max) < 1e-6 and abs(xmax - r_max) < 1e-6:
            edge_far.append(tag)
        elif abs(xmin) < 1e-6 and abs(xmax) < 1e-6 and tag not in line_filter:
            edge_axis.append(tag)

    # Physical Groups
    gmsh.model.addPhysicalGroup(2, surf_aq,    tag=1, name="aquifer")
    gmsh.model.addPhysicalGroup(2, surf_ct,    tag=2, name="caprock_top")
    gmsh.model.addPhysicalGroup(2, surf_cb,    tag=3, name="caprock_bottom")
    gmsh.model.addPhysicalGroup(1, line_filter, tag=4, name="filter_line")
    gmsh.model.addPhysicalGroup(0, [top_point], tag=5, name="filter_top")
    gmsh.model.addPhysicalGroup(1, edge_top,   tag=10, name="top")
    gmsh.model.addPhysicalGroup(1, edge_bot,   tag=11, name="bottom")
    gmsh.model.addPhysicalGroup(1, edge_far,   tag=12, name="far")

    # Mesh-Verfeinerung um die Filter-Linie
    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "CurvesList", line_filter)
    f_thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_thr, "InField", f_dist)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMin", m["size_in_filter_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMax", m["size_far_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "DistMin", m["well_size_radius_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "DistMax", m["well_size_radius_far_m"])
    gmsh.model.mesh.field.setAsBackgroundMesh(f_thr)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.model.mesh.generate(2)
    gmsh.write(str(msh_path))
    gmsh.finalize()
    return msh_path


# ======================================================================
#  Cycle-Curves (identisch zur Volumen-Variante)
# ======================================================================
def build_cycle_curves(cfg: dict) -> dict:
    cyc = cfg["cycles"]
    n = cyc["n_cycles"]; ramp = max(60.0, cyc["ramp_days"] * DAY)
    T0 = cfg["initial"]["T_K"]; T_hot = cfg["operation"]["T_hot_K"]
    rh = T0 / T_hot
    t_c  = cyc["charge_days"] * DAY
    t_sc = cyc["storage_after_charge_days"] * DAY
    t_d  = cyc["discharge_days"] * DAY
    t_sd = cyc["storage_after_discharge_days"] * DAY
    phases = [(t_c, +1.0, 1.0), (t_sc, 0.0, rh), (t_d, -1.0, rh), (t_sd, 0.0, rh)]
    times = [0.0]; v_m = [0.0]; v_T = [rh]; t_now = 0.0
    for _ in range(n):
        for dur, ms, ts in phases:
            if dur <= 0: continue
            t_now += ramp
            times.append(t_now); v_m.append(ms); v_T.append(ts)
            hold = max(0.0, dur - ramp)
            if hold > 0:
                t_now += hold; times.append(t_now); v_m.append(ms); v_T.append(ts)
    t_now += ramp; times.append(t_now); v_m.append(0.0); v_T.append(rh)
    return {"t_total": t_now,
            "cycle_mass": (np.array(times), np.array(v_m)),
            "cycle_T":    (np.array(times), np.array(v_T))}


# ======================================================================
#  PRJ
# ======================================================================
def _se(parent, tag, text=None, **attrs):
    el = ET.SubElement(parent, tag, **{k: str(v) for k, v in attrs.items()})
    if text is not None: el.text = str(text)
    return el
def _const_prop(parent, name, value):
    p = _se(parent, "property"); _se(p, "name", name); _se(p, "type", "Constant"); _se(p, "value", value)

def _add_phase_fluid(phases, fluid, op):
    ph = _se(phases, "phase"); _se(ph, "type", "AqueousLiquid")
    props = _se(ph, "properties")
    _const_prop(props, "density",                fluid["rho_ref_kg_m3"])
    _const_prop(props, "viscosity",              fluid["viscosity_Pa_s"])
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
    _add_phase_fluid(phases, fluid, op); _add_phase_solid(phases, mat, op)
    props = _se(med, "properties")
    _const_prop(props, "porosity",     mat["porosity"])
    _const_prop(props, "permeability", mat["permeability_m2"])
    p = _se(props, "property"); _se(p, "name", "thermal_conductivity")
    _se(p, "type", "EffectiveThermalConductivityPorosityMixing")
    _const_prop(props, "thermal_longitudinal_dispersivity", disp["alpha_L_m"])
    _const_prop(props, "thermal_transversal_dispersivity",  disp["alpha_T_m"])
    _const_prop(props, "storage", 0.0)

def _indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip(): elem.text = i + "  "
        for child in elem: _indent(child, level + 1)
        if not child.tail or not child.tail.strip(): child.tail = i
    if level and (not elem.tail or not elem.tail.strip()): elem.tail = i


def build_prj(cfg: dict, out_dir: Path, mesh_files: dict, curves: dict) -> Path:
    prefix = cfg["output"]["prefix"]
    fluid, op, init, sol, disp = (cfg["fluid"], cfg["operation"],
                                  cfg["initial"], cfg["solver"], cfg["dispersion"])
    n_steps = int(curves["t_total"] // cfg["time"]["dt_seconds"]) + 1

    root = ET.Element("OpenGeoSysProject")
    meshes = _se(root, "meshes")
    for key in ("domain", "top", "bottom", "far", "filter_line", "filter_top"):
        _se(meshes, "mesh", mesh_files[key], axially_symmetric="true")

    processes = _se(root, "processes")
    proc = _se(processes, "process")
    _se(proc, "name", "HT"); _se(proc, "type", "HT"); _se(proc, "integration_order", 2)
    pv = _se(proc, "process_variables"); _se(pv, "temperature", "T"); _se(pv, "pressure", "p")
    sv = _se(proc, "secondary_variables")
    _se(sv, "secondary_variable", internal_name="darcy_velocity", output_name="darcy_velocity")
    _se(proc, "specific_body_force", "0 0")

    # Media: aquifer = 0, caprock_top = 1, caprock_bottom = 2
    media = _se(root, "media")
    _add_medium(media, 0, cfg["materials"]["aquifer"],        fluid, op, disp)
    _add_medium(media, 1, cfg["materials"]["caprock_top"],    fluid, op, disp)
    _add_medium(media, 2, cfg["materials"]["caprock_bottom"], fluid, op, disp)

    # Time loop / output / solver — wie in der Volumen-Variante
    tl = _se(root, "time_loop")
    procs = _se(tl, "processes")
    pref = _se(procs, "process", ref="HT")
    nls = _se(pref, "nonlinear_solver"); _se(nls, "name", "basic_picard")
    convs = _se(pref, "convergence_criteria")
    for var, rtol in (("T", sol["rel_tol_T"]), ("p", sol["rel_tol_p"])):
        c = _se(convs, "convergence_criterion")
        _se(c, "type", "DeltaX"); _se(c, "norm_type", "NORM2"); _se(c, "reltol", rtol)
    ts = _se(pref, "time_stepping")
    _se(ts, "type", "FixedTimeStepping")
    _se(ts, "t_initial", 0); _se(ts, "t_end", curves["t_total"])
    tsteps = _se(ts, "timesteps"); pair = _se(tsteps, "pair")
    _se(pair, "repeat", n_steps); _se(pair, "delta_t", cfg["time"]["dt_seconds"])
    out_el = _se(tl, "output")
    _se(out_el, "type", "VTK"); _se(out_el, "prefix", prefix)
    out_steps = _se(out_el, "timesteps"); pair = _se(out_steps, "pair")
    _se(pair, "repeat", n_steps); _se(pair, "each_steps", cfg["time"]["output_every_n_steps"])

    # Parameters
    params = _se(root, "parameters")
    _se(_se(params, "parameter"), "name", "T0")  # filled below
    last = list(params)[-1]; _se(last, "type", "Constant"); _se(last, "value", init["T_K"])
    _se(_se(params, "parameter"), "name", "p0"); last = list(params)[-1]
    _se(last, "type", "Constant"); _se(last, "value", init["p_Pa"])
    # Massen-/T-Steuerung am Top-Knoten
    _se(_se(params, "parameter"), "name", "m_amp"); last = list(params)[-1]
    _se(last, "type", "Constant"); _se(last, "value", op["mass_flow_rate_kg_s"])
    _se(_se(params, "parameter"), "name", "T_hot"); last = list(params)[-1]
    _se(last, "type", "Constant"); _se(last, "value", op["T_hot_K"])
    # Curve-skaliert: NodalSourceTerm-Wert = cycle_mass * m_amp
    pp = _se(params, "parameter")
    _se(pp, "name", "m_nodal"); _se(pp, "type", "CurveScaled")
    _se(pp, "curve", "cycle_mass"); _se(pp, "parameter", "m_amp")
    pp = _se(params, "parameter")
    _se(pp, "name", "T_nodal"); _se(pp, "type", "CurveScaled")
    _se(pp, "curve", "cycle_T"); _se(pp, "parameter", "T_hot")

    cv = _se(root, "curves")
    for name, key in (("cycle_mass", "cycle_mass"), ("cycle_T", "cycle_T")):
        c = _se(cv, "curve"); _se(c, "name", name)
        t, v = curves[key]
        _se(c, "coords", " ".join(f"{x:.6e}" for x in t))
        _se(c, "values", " ".join(f"{x:.6e}" for x in v))

    # Process-Variables + BCs + Source-Terms
    pvars = _se(root, "process_variables")
    # Temperatur
    pv_T = _se(pvars, "process_variable"); _se(pv_T, "name", "T")
    _se(pv_T, "components", 1); _se(pv_T, "order", 1)
    _se(pv_T, "initial_condition", "T0")
    bcs_T = _se(pv_T, "boundary_conditions")
    for face in ("top", "bottom"):
        bc = _se(bcs_T, "boundary_condition")
        _se(bc, "mesh", Path(mesh_files[face]).stem)
        _se(bc, "type", "Dirichlet"); _se(bc, "parameter", "T0")
    # Dirichlet-T am Filter-Top-Punkt = T_inj-Curve
    bc = _se(bcs_T, "boundary_condition")
    _se(bc, "mesh", Path(mesh_files["filter_top"]).stem)
    _se(bc, "type", "Dirichlet"); _se(bc, "parameter", "T_nodal")

    # Druck
    pv_p = _se(pvars, "process_variable"); _se(pv_p, "name", "p")
    _se(pv_p, "components", 1); _se(pv_p, "order", 1)
    _se(pv_p, "initial_condition", "p0")
    bcs_p = _se(pv_p, "boundary_conditions")
    for face in ("top", "bottom", "far"):
        bc = _se(bcs_p, "boundary_condition")
        _se(bc, "mesh", Path(mesh_files[face]).stem)
        _se(bc, "type", "Dirichlet"); _se(bc, "parameter", "p0")
    # NodalSourceTerm am Filter-Top für Massen­injektion
    sts_p = _se(pv_p, "source_terms")
    st = _se(sts_p, "source_term")
    _se(st, "mesh", Path(mesh_files["filter_top"]).stem)
    _se(st, "type", "Nodal"); _se(st, "parameter", "m_nodal")

    # Solver
    nls = _se(root, "nonlinear_solvers")
    n = _se(nls, "nonlinear_solver")
    _se(n, "name", "basic_picard"); _se(n, "type", "Picard")
    _se(n, "max_iter", sol["nonlinear_iter"])
    _se(n, "linear_solver", "general_linear_solver")
    lss = _se(root, "linear_solvers")
    ls = _se(lss, "linear_solver"); _se(ls, "name", "general_linear_solver")
    eig = _se(ls, "eigen"); _se(eig, "solver_type", "BiCGSTAB")
    _se(eig, "precon_type", "ILUT"); _se(eig, "max_iteration_step", sol["linear_iter"])
    _se(eig, "error_tolerance", sol["linear_tol"]); _se(eig, "scaling", "true")

    _indent(root)
    prj_path = out_dir / f"{prefix}.prj"
    ET.ElementTree(root).write(prj_path, encoding="ISO-8859-1", xml_declaration=True)
    return prj_path


def run_ogs(prj_path: Path) -> int:
    ogs = shutil.which("ogs") or shutil.which("ogs.exe")
    if not ogs:
        print("ogs.exe nicht im PATH — Setup ist fertig, aber kein Lauf.")
        return 1
    return subprocess.call([ogs, str(prj_path), "-o", str(prj_path.parent)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-mesh", action="store_true")
    ap.add_argument("--no-run",  action="store_true")
    args = ap.parse_args()

    out_dir = Path(CONFIG["output"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = CONFIG["output"]["prefix"]
    msh_path = out_dir / f"{prefix}.msh"

    if not args.no_mesh:
        print("[1/3] gmsh: 2D-Radial-Mesh mit 1D-Filter-Linie ...")
        build_mesh(CONFIG, out_dir)
        print(f"      -> {msh_path}")
        print("[2/3] msh2vtu ...")
        msh2vtu(filename=msh_path, output_path=out_dir,
                output_prefix=prefix, dim=2, reindex=True)

    mesh_files = {
        "domain":      f"{prefix}_domain.vtu",
        "top":         f"{prefix}_physical_group_top.vtu",
        "bottom":      f"{prefix}_physical_group_bottom.vtu",
        "far":         f"{prefix}_physical_group_far.vtu",
        "filter_line": f"{prefix}_physical_group_filter_line.vtu",
        "filter_top":  f"{prefix}_physical_group_filter_top.vtu",
    }

    print("[3/3] OGS-Projektdatei (Linienquelle, Top-Knoten-Injektion) ...")
    curves = build_cycle_curves(CONFIG)
    prj_path = build_prj(CONFIG, out_dir, mesh_files, curves)
    print(f"      -> {prj_path}")

    if args.no_run: return 0
    return run_ogs(prj_path)


if __name__ == "__main__":
    sys.exit(main())
