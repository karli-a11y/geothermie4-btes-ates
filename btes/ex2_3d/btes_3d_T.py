#!/usr/bin/env python3
"""
BTES 3D — reine Waermeleitung (HEAT_CONDUCTION) auf dem Viertelmodell.

Rechnet dieselbe Aufgabe wie btes_3d.py, aber mit zwei Vereinfachungen, die
beide keine Naeherung sind, sondern aus der Aufgabenstellung selbst folgen:

    1) Prozess T statt HT      — die Druckgleichung ist hier trivial
    2) Viertel statt Vollfeld  — das Problem ist spiegelsymmetrisch

Zusammen sinkt die Zahl der Unbekannten auf ein Achtel. Beide Vereinfachungen
sind unten begruendet, samt der Bedingungen, unter denen sie NICHT mehr gelten.

-----------------------------------------------------------------------------
1) WARUM T STATT HT
-----------------------------------------------------------------------------
Der HT-Prozess loest Temperatur UND Druck. Im BTES-Aufbau ist das Druckfeld
aber identisch null, und zwar nicht naeherungsweise, sondern exakt:

  * kein thermischer Auftrieb   (fluid.beta_1_per_K = 0)
  * keine Schwerkraft           (specific_body_force = 0 0 0)
  * p = 0 als Dirichlet-Rand oben und unten, Anfangswert p = 0
  * an den Seitenflaechen keine Randbedingung

Damit ist die Druckgleichung homogen und von der Temperatur entkoppelt; ihre
einzige Loesung ist p = 0 ueberall und zu jeder Zeit. Aus p = 0 folgt eine
verschwindende Darcy-Geschwindigkeit, damit entfaellt der Advektionsterm der
Energiegleichung, und uebrig bleibt reine Waermeleitung. Wer den HT-Lauf
mitprotokolliert, sieht das direkt: das Konvergenzkriterium meldet fuer die
Druckkomponente |dx| = 0 und |x| = 0.

Auch wenn man einen regionalen Grundwassergradienten ansetzen WUERDE, bliebe
die Waermeleitung dominant. Mit einem hydraulischen Gradienten i = 2e-3 und der
durchlaessigsten Schicht dieses Schichtstapels (k ~ 5e-14 m2) ergibt sich

    v = k/mu * rho*g*i ~ 1e-9 m/s  ~  3 cm im Jahr

Die zugehoerige Peclet-Zahl ueber den Sondenabstand L = 6 m,

    Pe = v * rho_f * cp_f * L / lambda  ~  0.03

liegt zwei Groessenordnungen unter 1. Damit Advektion ueberhaupt mitspielt,
braeuchte es rund k ~ 4e-12 m2, also etwa das Achtzigfache — einen Sand- oder
Kiesaquifer, nicht die dichte Abfolge dieses Modells. Als Faustregel wird
Grundwasserstroemung fuer Erdwaermesonden ab etwa 1e-8 m/s relevant.

Die Stoffwerte muessen dafuer einmal vorab gemischt werden. HEAT_CONDUCTION
kennt keine Phasen, erwartet also je Medium direkt

    thermal_conductivity   = (1-n) * lambda_s + n * lambda_f
    density * cp           = (1-n) * rho_s*cp_s + n * rho_f*cp_f

Genau diese beiden Mischungen bildet der HT-Prozess intern ebenfalls
(EffectiveThermalConductivityPorosityMixing bzw. die volumetrische
Waermekapazitaet). Siehe effective_material(). In die Gleichung geht nur das
Produkt density*cp ein; die Aufteilung auf beide Groessen ist frei und hier so
gewaehlt, dass beide Zahlen physikalisch lesbar bleiben.

NICHT mehr anwendbar, sobald eine der Voraussetzungen faellt: eine
Druckrandbedingung an den Seitenflaechen (regionale Grundwasserstroemung), ein
wirklich durchlaessiger Aquifer im Schichtstapel, thermischer Auftrieb oder
Schwerkraft. Dann gehoert der HT-Prozess zurueck ins Modell.

-----------------------------------------------------------------------------
2) WARUM DAS VIERTELMODELL GENUEGT
-----------------------------------------------------------------------------
Das Sondenfeld ist ein regelmaessiges Rechteckraster, mittig im Gebiet. Alle
Sonden tragen dieselbe Lastkurve, die Schichten liegen waagerecht, und die
Randbedingungen oben und unten sind gleichfoermig. Damit sind die beiden
senkrechten Ebenen x = 0 und y = 0 Spiegelebenen des GESAMTEN Problems —
Geometrie, Material, Last und Randbedingungen zugleich.

Auf einer Spiegelebene verschwindet der Waermestrom senkrecht zur Ebene. In OGS
ist eine Flaeche OHNE Randbedingung genau das: ein Nullfluss-Rand. Die
Schnittflaechen bekommen deshalb bewusst keine Randbedingung — damit ist die
Symmetriebedingung exakt erfuellt, ohne dass etwas hinzugefuegt werden muesste.
Das Viertel liefert dieselbe Loesung wie das Vollfeld, nicht eine genaeherte.

Wichtig ist dabei, dass gerade der interessante Teil des Speichers im Modell
bleibt. Die staerkste Abkuehlung bzw. Aufheizung tritt in der Feldmitte auf,
weil sich dort die Nachbarsonden am staerksten gegenseitig beeinflussen; die
Randsonden sind stets die unkritischen. Die Feldmitte liegt aber genau auf dem
Schnittpunkt der beiden Symmetrieebenen und gehoert damit zum Viertel. Weggelassen
wird nur, was sich ohnehin spiegelbildlich wiederholt.

Voraussetzung ist eine GERADE Sondenzahl je Richtung. Dann liegt keine Sonde auf
einer Schnittebene und jede Sondenbox bleibt vollstaendig im Modell. Bei
ungerader Zahl saesse eine Sonde halbiert auf der Ebene und ihr Quellterm
muesste mit halbiert werden — das Skript bricht in diesem Fall ab, statt still
etwas Falsches zu rechnen.

Die AEUSSEREN Seitenflaechen behalten die Behandlung des Vollmodells (hier:
keine Randbedingung, also adiabat). Die Viertelung aendert daran nichts.

NICHT anwendbar bei unsymmetrischem Feldgrundriss, bei einzeln geregelten
Sonden mit unterschiedlichen Lasten, oder bei regionaler Grundwasserstroemung —
letztere bricht die Symmetrie in Stroemungsrichtung, ueblich ist dann noch ein
Halbmodell mit Schnitt quer zur Stroemung.

-----------------------------------------------------------------------------
3) WIE DER WAERMEBEDARF ANZUGEBEN IST
-----------------------------------------------------------------------------
Angegeben wird die GESAMTLAST DES FELDES, nicht die Last je Sonde. Das Skript
teilt intern durch die Sondenzahl des VOLLEN Feldes:

    P_je_Sonde = field_power_kW * 1000 / (n_x_full * n_y_full)

und praegt diesen Wert jeder im Viertel modellierten Sonde auf. Jede Sonde
traegt damit exakt dieselbe Leistung wie im Vollmodell.

Das ist die haeufigste Fehlerquelle beim Viertelmodell: wer die Feldlast durch
die Sondenzahl des VIERTELS teilt, rechnet mit der vierfachen Last. Deshalb
nimmt dieses Skript die Feldlast entgegen und teilt selbst.

Vorzeichen: positiv = Beladung (Waerme in den Untergrund), negativ = Foerderung.
Der Quellterm wirkt volumetrisch auf der Sondenbox,

    q_v = P_je_Sonde / (borehole_dx_m * borehole_dy_m * Sondenlaenge)

Monatsprofil: zwoelf Werte in kW, Feld gesamt. Beispiel — ein Feld, das im
Januar 1.200 kW liefern soll und im Juli 700 kW einspeichert, bekommt bei
220 Sonden je Sonde -5.455 W bzw. +3.182 W. Als Kontrolle druckt das Skript die
spezifische Entzugsrate in W/m; der Literaturbereich fuer Erdwaermesonden liegt
bei 20 bis 70 W/m.

-----------------------------------------------------------------------------
Aufruf:
    python btes_3d_T.py            # Mesh + .prj + OGS-Lauf
    python btes_3d_T.py --no-run   # nur Setup, kein OGS
    python btes_3d_T.py --no-mesh  # nur .prj erzeugen
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


def _safe_name(name):
    """Physical-Group-Namen dateisystemsicher machen (wie in btes_3d.py)."""
    keep = ("abcdefghijklmnopqrstuvwxyz"
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")
    return "".join(c if c in keep else "_" for c in str(name))


# --- ogstools >=0.8 Kompatibilitaets-Shim fuer die alte msh2vtu-API ---
def msh2vtu(filename, output_path, output_prefix, dim, reindex=True):
    import ogstools as ot
    from pathlib import Path as _P
    meshes = ot.Meshes.from_gmsh(filename=str(filename), dim=dim,
                                 reindex=reindex, log=False)
    output_path = _P(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    for name, mesh in meshes.items():
        fname = (f"{output_prefix}_domain.vtu" if name == "domain"
                 else f"{output_prefix}_physical_group_{_safe_name(name)}.vtu")
        mesh.save(str(output_path / fname), binary=True)


# ======================================================================
#  CONFIG  --  hier alles anpassen
# ======================================================================
CONFIG: dict = {
    "domain": {
        # Im Viertelmodell sind das die Abmessungen des QUADRANTEN, das
        # entsprechende Vollgebiet ist doppelt so gross. Bei symmetry="full"
        # bezeichnen sie die halbe Kantenlaenge, das Gebiet ist dann
        # [-size_x_m, +size_x_m] x [-size_y_m, +size_y_m] — in beiden Faellen
        # beschreibt dieselbe Zahl also dasselbe Vollgebiet.
        "size_x_m":  60.0,
        "size_y_m":  60.0,
        "z_base_m":   0.0,
    },
    # ------------------------------------------------------------------
    # SCHICHTEN (von OBEN nach UNTEN), wie in btes_3d.py
    # ------------------------------------------------------------------
    "layers": [
        {"name": "Deckschicht", "thickness_m": 15.0,
         "porosity": 0.15, "rho_s_kg_m3": 2500.0,
         "cp_s_J_kgK": 880.0, "lambda_s_W_mK": 2.2},
        {"name": "Sandstein", "thickness_m": 60.0,
         "porosity": 0.20, "rho_s_kg_m3": 2650.0,
         "cp_s_J_kgK": 560.0, "lambda_s_W_mK": 2.5},
        {"name": "Grundgebirge", "thickness_m": 45.0,
         "porosity": 0.05, "rho_s_kg_m3": 2700.0,
         "cp_s_J_kgK": 750.0, "lambda_s_W_mK": 2.8},
    ],
    "borehole": {
        "depth_top_m":     2.0,
        "depth_bottom_m": 100.0,
        "borehole_dx_m":   0.6,
        "borehole_dy_m":   0.6,
        # nur fuer die Nachrechnung der Fluidtemperatur, nicht fuer die Loesung
        "r_borehole_m":    0.075,
        "R_b_Km_per_W":    0.10,
    },
    "field": {
        # Sondenzahl des VOLLEN Feldes. Im Viertelmodell muessen beide
        # Zahlen gerade sein, damit keine Sonde auf einer Schnittebene liegt.
        "n_x_full":   8,
        "n_y_full":   8,
        "spacing_m":  6.0,
        "symmetry":   "quarter",     # "quarter" oder "full"
    },
    "mesh": {
        "size_in_borehole_m":      0.4,
        "size_near_field_m":       1.5,
        "size_far_m":             12.0,
        "field_size_radius_m":     8.0,
        "field_size_radius_far_m": 30.0,
    },
    "fluid": {
        # Porenwasser. Geht nur ueber die Mischung in die Stoffwerte ein —
        # ohne Stroemung gibt es keinen Transport.
        "rho_ref_kg_m3": 1000.0,
        "cp_J_kgK":      4180.0,
        "lambda_W_mK":      0.6,
    },
    "initial": {
        "T_K":                         283.15,
        "T_surface_K":                 283.15,
        "geothermal_gradient_K_per_m":   0.0,   # 0.03 fuer realistischen Gradient
    },
    "operation": {
        # Referenzleistung des GESAMTEN Feldes [kW]. Dient als Bezugsgroesse,
        # auf die die Lastkurve skaliert wird. Aufteilung auf die Sonden
        # geschieht intern, siehe Kopf, Abschnitt 3.
        "field_power_kW": 400.0,
    },
    # ------------------------------------------------------------------
    # ZYKLEN
    # ------------------------------------------------------------------
    # Modus A: 4-Phasen-Zyklus (Beladung, Pause, Foerderung, Pause).
    # Modus B: cycles.monthly_field_power_kW auf zwoelf Monatswerte setzen,
    #          Feld GESAMT in kW, positiv = Beladung. Ueberschreibt Modus A.
    # ------------------------------------------------------------------
    "cycles": {
        "n_cycles":                     1,
        "charge_days":                 91.25,
        "storage_after_charge_days":   91.25,
        "discharge_days":              91.25,
        "storage_after_discharge_days": 91.25,
        "ramp_days":                    7.0,
        "monthly_field_power_kW":      None,
    },
    "time": {
        "dt_seconds":           7 * 86400.0,
        "output_every_n_steps": 1,
    },
    "output": {
        "prefix":    "btes_3d_T",
        "out_dir":   "out",
        "variables": ["T"],
    },
    "solver": {
        # Reine Waermeleitung -> Systemmatrix symmetrisch positiv definit,
        # also CG statt BiCGSTAB. ILUT waere hier unnoetig teuer.
        "solver_type":    "CG",
        "precon_type":    "DIAGONAL",
        "linear_tol":     1.0e-10,
        "linear_iter":    20000,
        "nonlinear_iter": 20,
        "rel_tol_T":      1.0e-2,
    },
}

DAY = 86400.0


# ======================================================================
#  Geometrie
# ======================================================================
def _is_quarter(cfg: dict) -> bool:
    return cfg["field"].get("symmetry", "quarter") == "quarter"


def _layer_stack(cfg: dict):
    """Schichtliste von oben nach unten in z-Grenzen umrechnen (bottom-up)."""
    z = cfg["domain"].get("z_base_m", 0.0)
    out = []
    for L in reversed(list(cfg["layers"])):
        out.append({**L, "z_low": z, "z_high": z + float(L["thickness_m"])})
        z += float(L["thickness_m"])
    return out, z


def _n_boreholes_full(cfg: dict) -> int:
    return int(cfg["field"]["n_x_full"]) * int(cfg["field"]["n_y_full"])


def _borehole_positions(cfg: dict) -> list[tuple[float, float]]:
    """Sondenpositionen; im Viertelmodell nur der Quadrant x > 0, y > 0."""
    nx, ny = int(cfg["field"]["n_x_full"]), int(cfg["field"]["n_y_full"])
    s = float(cfg["field"]["spacing_m"])
    xs = [-(nx - 1) * s / 2.0 + i * s for i in range(nx)]
    ys = [-(ny - 1) * s / 2.0 + i * s for i in range(ny)]
    if not _is_quarter(cfg):
        return [(x, y) for x in xs for y in ys]
    if nx % 2 or ny % 2:
        raise ValueError(
            f"Viertelmodell verlangt gerade Sondenzahlen je Richtung "
            f"(n_x_full={nx}, n_y_full={ny}). Bei ungerader Zahl laege eine "
            f"Sonde auf einer Symmetrieebene und muesste halbiert werden. "
            f"Entweder Sondenzahl aendern oder symmetry='full' setzen."
        )
    return [(x, y) for x in xs if x > 1e-9 for y in ys if y > 1e-9]


def _power_per_borehole_W(cfg: dict) -> float:
    """Referenzleistung je Sonde [W] aus der Feldgesamtlast."""
    return cfg["operation"]["field_power_kW"] * 1000.0 / _n_boreholes_full(cfg)


# ======================================================================
#  Effektive Stoffwerte  (Korn + Porenwasser -> Ersatzmedium)
# ======================================================================
def effective_material(mat: dict, fluid: dict) -> dict:
    """Mischt so, wie der HT-Prozess es intern ebenfalls tut.

    lambda_eff = (1-n)*lambda_s + n*lambda_f   (PorosityMixing)
    (rho*cp)_eff = (1-n)*rho_s*cp_s + n*rho_f*cp_f

    In die Waermeleitungsgleichung geht nur das Produkt density*cp ein. Die
    Dichte wird deshalb als Volumenmittel angesetzt und cp so bestimmt, dass
    das Produkt stimmt — beide Zahlen bleiben damit physikalisch lesbar.
    """
    n = float(mat["porosity"])
    rho_s, cp_s = mat["rho_s_kg_m3"], mat["cp_s_J_kgK"]
    rho_f, cp_f = fluid["rho_ref_kg_m3"], fluid["cp_J_kgK"]
    rho_eff = (1.0 - n) * rho_s + n * rho_f
    rhoc_eff = (1.0 - n) * rho_s * cp_s + n * rho_f * cp_f
    return {
        "density":                rho_eff,
        "specific_heat_capacity": rhoc_eff / rho_eff,
        "thermal_conductivity":   (1.0 - n) * mat["lambda_s_W_mK"]
                                  + n * fluid["lambda_W_mK"],
        "rho_c":                  rhoc_eff,
    }


# ======================================================================
#  Netz — gmsh
# ======================================================================
def build_mesh(cfg: dict, out_dir: Path) -> Path:
    msh_path = out_dir / f"{cfg['output']['prefix']}.msh"

    LX, LY = cfg["domain"]["size_x_m"], cfg["domain"]["size_y_m"]
    z_base = cfg["domain"]["z_base_m"]
    layers, z_top = _layer_stack(cfg)

    quarter = _is_quarter(cfg)
    x_lo, y_lo = (0.0, 0.0) if quarter else (-LX, -LY)
    x_hi, y_hi = LX, LY
    Wx, Wy = x_hi - x_lo, y_hi - y_lo

    bh = cfg["borehole"]
    dx_b, dy_b = bh["borehole_dx_m"], bh["borehole_dy_m"]
    z_bh_top = z_top - bh["depth_top_m"]
    z_bh_bot = z_top - bh["depth_bottom_m"]
    if z_bh_top <= z_bh_bot:
        raise ValueError("depth_bottom_m muss groesser als depth_top_m sein.")
    if z_bh_bot < z_base - 1e-9 or z_bh_top > z_top + 1e-9:
        raise ValueError("Sondentiefe liegt ausserhalb der Schichtdomaene.")

    pos = _borehole_positions(cfg)
    for x, y in pos:
        if not (x_lo < x - dx_b / 2 and x + dx_b / 2 < x_hi
                and y_lo < y - dy_b / 2 and y + dy_b / 2 < y_hi):
            raise ValueError(
                f"Sondenbox bei ({x:.1f}, {y:.1f}) passt nicht ins Gebiet — "
                f"domain.size_x_m / size_y_m vergroessern.")

    mh = cfg["mesh"]
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("btes_T")

    layer_boxes = [gmsh.model.occ.addBox(x_lo, y_lo, L["z_low"], Wx, Wy,
                                         L["z_high"] - L["z_low"])
                   for L in layers]
    bh_boxes = [gmsh.model.occ.addBox(x - dx_b / 2, y - dy_b / 2, z_bh_bot,
                                      dx_b, dy_b, z_bh_top - z_bh_bot)
                for x, y in pos]
    gmsh.model.occ.fragment([(3, layer_boxes[0])],
                            [(3, b) for b in layer_boxes[1:]]
                            + [(3, b) for b in bh_boxes])
    gmsh.model.occ.synchronize()

    # Volumen zuordnen: kleine Grundflaeche im Sonden-z-Bereich -> Sonde
    bh_pos = np.array(pos)
    vol_bh = {i: [] for i in range(len(pos))}
    vol_layer = {i: [] for i in range(len(layers))}
    for _, tag in gmsh.model.getEntities(3):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(3, tag)
        zc = 0.5 * (zmin + zmax)
        xc, yc = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)
        if (xmax - xmin) < 0.3 * Wx and (z_bh_bot - 1e-3) <= zc <= (z_bh_top + 1e-3):
            d = np.hypot(bh_pos[:, 0] - xc, bh_pos[:, 1] - yc)
            vol_bh[int(np.argmin(d))].append(tag)
            continue
        for i, L in enumerate(layers):
            if L["z_low"] - 1e-6 <= zc <= L["z_high"] + 1e-6:
                vol_layer[i].append(tag)
                break
    missing = [i for i, t in vol_bh.items() if not t]
    if missing:
        raise RuntimeError(f"Sonden ohne Volumen gefunden (Index {missing}).")

    # Aussenflaechen; im Viertelmodell zusaetzlich die Symmetrieebenen
    def _flat(lo, hi, v):
        return abs(lo - v) < 1e-6 and abs(hi - v) < 1e-6

    surf_top, surf_bot, surf_out, surf_sym = [], [], [], []
    for _, tag in gmsh.model.getEntities(2):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(2, tag)
        zc = 0.5 * (zmin + zmax)
        if (xmax - xmin) >= 0.9 * Wx and abs(zc - z_top) < 1e-6:
            surf_top.append(tag); continue
        if (xmax - xmin) >= 0.9 * Wx and abs(zc - z_base) < 1e-6:
            surf_bot.append(tag); continue
        on_sym = quarter and (_flat(xmin, xmax, 0.0) or _flat(ymin, ymax, 0.0))
        on_out = (_flat(xmin, xmax, x_hi) or _flat(ymin, ymax, y_hi)
                  or (not quarter and (_flat(xmin, xmax, x_lo)
                                       or _flat(ymin, ymax, y_lo))))
        if on_sym:
            surf_sym.append(tag)
        elif on_out:
            surf_out.append(tag)

    pg = 1
    for i, L in enumerate(layers):
        gmsh.model.addPhysicalGroup(3, vol_layer[i], tag=pg, name=L["name"]); pg += 1
    for i in range(len(pos)):
        gmsh.model.addPhysicalGroup(3, vol_bh[i], tag=pg, name=f"bh_{i:02d}"); pg += 1
    gmsh.model.addPhysicalGroup(2, surf_top, tag=100, name="top")
    gmsh.model.addPhysicalGroup(2, surf_bot, tag=101, name="bottom")
    if surf_out:
        gmsh.model.addPhysicalGroup(2, surf_out, tag=102, name="lateral")
    if surf_sym:
        # Wird bewusst NICHT als Randbedingung verwendet: eine Flaeche ohne
        # Randbedingung ist in OGS ein Nullfluss-Rand, und das ist genau die
        # Symmetriebedingung. Die Gruppe existiert nur zur Kontrolle.
        gmsh.model.addPhysicalGroup(2, surf_sym, tag=103, name="symmetry")

    all_bh = [t for tags in vol_bh.values() for t in tags]
    bh_surf = []
    for tag in all_bh:
        for d, t in gmsh.model.getBoundary([(3, tag)], oriented=False):
            if d == 2:
                bh_surf.append(abs(t))
    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "SurfacesList", sorted(set(bh_surf)))
    f_thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_thr, "InField", f_dist)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMin", mh["size_near_field_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMax", mh["size_far_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "DistMin", mh["field_size_radius_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "DistMax", mh["field_size_radius_far_m"])
    gmsh.model.mesh.field.setAsBackgroundMesh(f_thr)

    bh_pts = []
    for tag in all_bh:
        for d, t in gmsh.model.getBoundary([(3, tag)], recursive=True,
                                           oriented=False):
            if d == 0:
                bh_pts.append((d, t))
    if bh_pts:
        gmsh.model.mesh.setSize(bh_pts, mh["size_in_borehole_m"])

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.model.mesh.generate(3)
    gmsh.write(str(msh_path))
    n_nodes = len(gmsh.model.mesh.getNodes()[0])
    gmsh.finalize()
    print(f"      {len(pos)} Sonden im Modell, {n_nodes:,} Knoten, "
          f"{n_nodes:,} Freiheitsgrade (nur T)")
    return msh_path


def _mesh_files(cfg: dict) -> dict:
    p = cfg["output"]["prefix"]
    files = {"domain": f"{p}_domain.vtu",
             "top": f"{p}_physical_group_top.vtu",
             "bottom": f"{p}_physical_group_bottom.vtu",
             "_n_bh": len(_borehole_positions(cfg))}
    for L in cfg["layers"]:
        files[L["name"]] = f"{p}_physical_group_{_safe_name(L['name'])}.vtu"
    for i in range(files["_n_bh"]):
        files[f"bh_{i:02d}"] = f"{p}_physical_group_bh_{i:02d}.vtu"
    return files


def convert_mesh(cfg: dict, msh_path: Path, out_dir: Path) -> dict:
    msh2vtu(msh_path, out_dir, cfg["output"]["prefix"], dim=3, reindex=True)
    return _mesh_files(cfg)


# ======================================================================
#  Lastkurve
# ======================================================================
def build_cycle_curves(cfg: dict) -> dict:
    """Eine Kurve fuer alle Sonden: +1 Beladung, 0 Pause, -1 Foerderung.

    Modus B skaliert die Monatswerte auf die Referenzleistung je Sonde. Weil
    beide aus derselben Feldgesamtlast abgeleitet werden, ist das Verhaeltnis
    unabhaengig von der Sondenzahl.
    """
    cyc = cfg["cycles"]
    ramp = max(60.0, cyc["ramp_days"] * DAY)
    n = cyc["n_cycles"]

    monthly = cyc.get("monthly_field_power_kW")
    if monthly is not None:
        assert len(monthly) == 12, \
            "cycles.monthly_field_power_kW muss 12 Werte enthalten."
        P_ref_kW = cfg["operation"]["field_power_kW"]
        if P_ref_kW == 0:
            raise ValueError("operation.field_power_kW muss > 0 sein.")
        month = 365.25 / 12.0 * DAY
        times, vals, t = [0.0], [0.0], 0.0
        for _ in range(n):
            for P in monthly:
                q = float(P) / P_ref_kW
                t += ramp; times.append(t); vals.append(q)
                hold = max(0.0, month - ramp)
                if hold > 0.0:
                    t += hold; times.append(t); vals.append(q)
        t += ramp; times.append(t); vals.append(0.0)
        return {"t_total": t, "cycle_q": (np.array(times), np.array(vals))}

    phases = [(cyc["charge_days"], +1.0),
              (cyc["storage_after_charge_days"], 0.0),
              (cyc["discharge_days"], -1.0),
              (cyc["storage_after_discharge_days"], 0.0)]
    times, vals, t = [0.0], [0.0], 0.0
    for _ in range(n):
        for dur, q in phases:
            if dur <= 0.0:
                continue
            t += ramp; times.append(t); vals.append(q)
            hold = max(0.0, dur * DAY - ramp)
            if hold > 0.0:
                t += hold; times.append(t); vals.append(q)
    t += ramp; times.append(t); vals.append(0.0)
    return {"t_total": t, "cycle_q": (np.array(times), np.array(vals))}


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
    _se(p, "name", name); _se(p, "type", "Constant"); _se(p, "value", value)


def _prop(parent, name, value):
    p = _se(parent, "property")
    _se(p, "name", name); _se(p, "type", "Constant"); _se(p, "value", value)


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
    fluid, init, sol = cfg["fluid"], cfg["initial"], cfg["solver"]
    bh = cfg["borehole"]

    h_bh = bh["depth_bottom_m"] - bh["depth_top_m"]
    V_bh = bh["borehole_dx_m"] * bh["borehole_dy_m"] * h_bh
    q_v = _power_per_borehole_W(cfg) / V_bh

    root = ET.Element("OpenGeoSysProject")

    meshes = _se(root, "meshes")
    keys = ["domain", "top", "bottom"] + [L["name"] for L in cfg["layers"]]
    keys += [f"bh_{i:02d}" for i in range(mesh_files["_n_bh"])]
    for k in keys:
        if mesh_files.get(k):
            _se(meshes, "mesh", mesh_files[k])

    proc = _se(_se(root, "processes"), "process")
    _se(proc, "name", "HeatConduction")
    _se(proc, "type", "HEAT_CONDUCTION")
    _se(proc, "integration_order", 2)
    _se(_se(proc, "process_variables"), "process_variable", "T")

    # Media: erst Schichten (unten -> oben, passend zu build_mesh), dann Sonden
    media = _se(root, "media")
    layers_bu = list(reversed(cfg["layers"]))
    for mid, L in enumerate(layers_bu):
        eff = effective_material(L, fluid)
        med = _se(media, "medium", id=mid)
        _se(med, "phases")
        props = _se(med, "properties")
        _prop(props, "thermal_conductivity", eff["thermal_conductivity"])
        _prop(props, "density", eff["density"])
        _prop(props, "specific_heat_capacity", eff["specific_heat_capacity"])
    # Sondenmaterial: mittlere Schicht (didaktische Voreinstellung; fuer eine
    # realistische Rechnung gehoert hier ein Verfuellmaterial hinein)
    eff_bh = effective_material(layers_bu[len(layers_bu) // 2], fluid)
    for i in range(mesh_files["_n_bh"]):
        med = _se(media, "medium", id=len(layers_bu) + i)
        _se(med, "phases")
        props = _se(med, "properties")
        _prop(props, "thermal_conductivity", eff_bh["thermal_conductivity"])
        _prop(props, "density", eff_bh["density"])
        _prop(props, "specific_heat_capacity", eff_bh["specific_heat_capacity"])

    tl = _se(root, "time_loop")
    pref = _se(_se(tl, "processes"), "process", ref="HeatConduction")
    _se(pref, "nonlinear_solver", "basic_picard")
    cc = _se(pref, "convergence_criterion")
    _se(cc, "type", "DeltaX"); _se(cc, "norm_type", "NORM2")
    _se(cc, "reltol", sol["rel_tol_T"])
    _se(_se(pref, "time_discretization"), "type", "BackwardEuler")
    ts = _se(pref, "time_stepping")
    _se(ts, "type", "FixedTimeStepping")
    _se(ts, "t_initial", 0.0); _se(ts, "t_end", curves["t_total"])
    n_steps = int(np.ceil(curves["t_total"] / cfg["time"]["dt_seconds"]))
    pair = _se(_se(ts, "timesteps"), "pair")
    _se(pair, "repeat", n_steps); _se(pair, "delta_t", cfg["time"]["dt_seconds"])

    out = _se(tl, "output")
    _se(out, "type", "VTK"); _se(out, "prefix", prefix)
    pair = _se(_se(out, "timesteps"), "pair")
    _se(pair, "repeat", n_steps)
    _se(pair, "each_steps", cfg["time"]["output_every_n_steps"])
    _se(out, "output_iteration_results", "false")
    ve = _se(out, "variables")
    for v in cfg["output"]["variables"]:
        _se(ve, "variable", v)

    params = _se(root, "parameters")
    grad = init.get("geothermal_gradient_K_per_m", 0.0)
    if abs(grad) > 1e-12:
        z_tot = sum(L["thickness_m"] for L in cfg["layers"])
        p_el = _se(params, "parameter")
        _se(p_el, "name", "T0"); _se(p_el, "type", "Function")
        _se(p_el, "expression",
            f"{init.get('T_surface_K', init['T_K'])} + ({grad:.6g})*({z_tot} - z)")
    else:
        _const_param(params, "T0", init["T_K"])
    _const_param(params, "q_v_amp", q_v)
    p = _se(params, "parameter")
    _se(p, "name", "q_v_borehole"); _se(p, "type", "CurveScaled")
    _se(p, "curve", "cycle_q"); _se(p, "parameter", "q_v_amp")

    c = _se(_se(root, "curves"), "curve")
    t_, v_ = curves["cycle_q"]
    _se(c, "name", "cycle_q")
    _se(c, "coords", " ".join(f"{x:.6e}" for x in t_))
    _se(c, "values", " ".join(f"{x:.6e}" for x in v_))

    pv = _se(_se(root, "process_variables"), "process_variable")
    _se(pv, "name", "T"); _se(pv, "components", 1); _se(pv, "order", 1)
    _se(pv, "initial_condition", "T0")
    bcs = _se(pv, "boundary_conditions")
    for face in ("top", "bottom"):
        bc = _se(bcs, "boundary_condition")
        _se(bc, "mesh", Path(mesh_files[face]).stem)
        _se(bc, "type", "Dirichlet"); _se(bc, "parameter", "T0")
    # Seitenflaechen bleiben ohne Randbedingung. Auf den Symmetrieebenen ist
    # das exakt die Spiegelbedingung, aussen wirkt es adiabat.
    sts = _se(pv, "source_terms")
    for i in range(mesh_files["_n_bh"]):
        st = _se(sts, "source_term")
        _se(st, "mesh", Path(mesh_files[f"bh_{i:02d}"]).stem)
        _se(st, "type", "Volumetric"); _se(st, "parameter", "q_v_borehole")

    nl = _se(_se(root, "nonlinear_solvers"), "nonlinear_solver")
    _se(nl, "name", "basic_picard"); _se(nl, "type", "Picard")
    _se(nl, "max_iter", sol["nonlinear_iter"])
    _se(nl, "linear_solver", "general_linear_solver")
    ls = _se(_se(root, "linear_solvers"), "linear_solver")
    _se(ls, "name", "general_linear_solver")
    eig = _se(ls, "eigen")
    _se(eig, "solver_type", sol["solver_type"])
    _se(eig, "precon_type", sol["precon_type"])
    _se(eig, "max_iteration_step", sol["linear_iter"])
    _se(eig, "error_tolerance", sol["linear_tol"])
    _se(eig, "scaling", "true")

    _indent(root)
    prj = out_dir / f"{prefix}.prj"
    ET.ElementTree(root).write(prj, encoding="ISO-8859-1", xml_declaration=True)
    return prj


# ======================================================================
#  Run + CLI
# ======================================================================
def run_ogs(prj_path: Path) -> int:
    exe = shutil.which("ogs") or shutil.which("ogs.exe")
    if not exe:
        print("ogs.exe nicht im PATH", file=sys.stderr)
        return 1
    cmd = [exe, str(prj_path), "-o", str(prj_path.parent)]
    print(">>", " ".join(cmd))
    return subprocess.call(cmd)


def _report(cfg: dict, curves: dict) -> None:
    """Kontrollausgabe der Lastaufteilung — die haeufigste Fehlerquelle."""
    n_full = _n_boreholes_full(cfg)
    n_mod = len(_borehole_positions(cfg))
    P_bh = _power_per_borehole_W(cfg)
    bh = cfg["borehole"]
    L = bh["depth_bottom_m"] - bh["depth_top_m"]
    monthly = cfg["cycles"].get("monthly_field_power_kW")
    peak_kW = (max(abs(p) for p in monthly) if monthly
               else abs(cfg["operation"]["field_power_kW"]))
    print(f"      Feld gesamt {n_full} Sonden, davon {n_mod} im Modell "
          f"({'Viertel' if _is_quarter(cfg) else 'Vollfeld'})")
    print(f"      Referenzlast {cfg['operation']['field_power_kW']:,.0f} kW Feld "
          f"-> {P_bh:,.0f} W je Sonde")
    print(f"      Spitzenlast {peak_kW:,.0f} kW Feld -> "
          f"{peak_kW*1000/n_full/L:.1f} W/m spezifisch "
          f"(Literatur Erdwaermesonden 20-70 W/m)")
    print(f"      t_end = {curves['t_total']/DAY:.1f} d")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="BTES 3D, reine Waermeleitung, Viertelmodell")
    ap.add_argument("--no-mesh", action="store_true")
    ap.add_argument("--no-run", action="store_true")
    args = ap.parse_args()

    out_dir = Path(CONFIG["output"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = CONFIG["output"]["prefix"]

    if not args.no_mesh:
        print("[1/3] gmsh: Sondenfeld ...")
        build_mesh(CONFIG, out_dir)
        print("[2/3] msh2vtu: Konvertierung ...")
        mesh_files = convert_mesh(CONFIG, out_dir / f"{prefix}.msh", out_dir)
    else:
        mesh_files = _mesh_files(CONFIG)

    print("[3/3] OGS-Projektdatei ...")
    curves = build_cycle_curves(CONFIG)
    prj = build_prj(CONFIG, out_dir, mesh_files, curves)
    print(f"      {prj}")
    _report(CONFIG, curves)

    if args.no_run:
        return 0
    print(">>> OGS starten")
    return run_ogs(prj)


if __name__ == "__main__":
    sys.exit(main())
