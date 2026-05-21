#!/usr/bin/env python3
"""
Auswertung & Plots für die ATES-3D-Übung (Doublet HW/CW).

Erzeugt in figures/:
  - 1_well_temperature.png       T(t) an HW und CW (Aquifer-Mitte)
  - 2_field_snapshots.png        T-Feld in der Aquifer-Mittel­ebene zu 4 Zeitpunkten
  - 3_energy_balance.png         Gespeicherte Wärmemenge & Recovery-Effizienz
  - 4_plume_extent.png           Reichweite der thermischen Front um HW
"""
from __future__ import annotations

import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

from pathlib import Path
from xml.etree import ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

from ates_3d import CONFIG

DELTA_T_THRESHOLD = 1.0
N_SNAPSHOTS       = 4
DAY = 86400.0


def read_pvd(pvd_path: Path):
    tree = ET.parse(pvd_path); items = []
    base = pvd_path.parent
    for ds in tree.iter("DataSet"):
        items.append((float(ds.attrib["timestep"]), base / ds.attrib["file"]))
    items.sort(); return items


def effective_rho_cp(cfg):
    aq = cfg["materials"]["aquifer"]
    phi = aq["porosity"]
    return (phi * cfg["fluid"]["rho_ref_kg_m3"] * cfg["fluid"]["cp_J_kgK"]
            + (1 - phi) * aq["rho_s_kg_m3"] * aq["cp_s_J_kgK"])


def cartesian_integral(mesh, field):
    vols = mesh.compute_cell_sizes(length=False, area=False, volume=True)["Volume"]
    if field.shape[0] == mesh.n_points:
        cell_field = np.empty(mesh.n_cells)
        for i in range(mesh.n_cells):
            pts = mesh.get_cell(i).point_ids
            cell_field[i] = field[pts].mean()
    else:
        cell_field = field
    return float(np.sum(cell_field * vols))


def main() -> int:
    out_dir = Path(CONFIG["output"]["out_dir"])
    prefix  = CONFIG["output"]["prefix"]
    pvd     = out_dir / f"{prefix}.pvd"
    if not pvd.exists():
        print(f"FEHLER: {pvd} nicht gefunden. Bitte zuerst die Simulation ausführen.")
        return 1

    fig_dir = Path(__file__).parent / "figures"
    fig_dir.mkdir(exist_ok=True)

    T0 = CONFIG["initial"]["T_K"]
    rho_cp = effective_rho_cp(CONFIG)

    hw_xy = CONFIG["wells"]["hot_well_xy"]
    cw_xy = CONFIG["wells"]["cold_well_xy"]
    z_aq_bot = CONFIG["layers"]["caprock_bottom_thickness_m"]
    z_aq_top = z_aq_bot + CONFIG["layers"]["aquifer_thickness_m"]
    z_mid    = 0.5 * (z_aq_bot + z_aq_top)

    steps = read_pvd(pvd)
    times = np.array([t for t, _ in steps]); times_d = times / DAY

    # 1) T(t) an HW/CW -------------------------------------------------
    probes = np.array([[hw_xy[0], hw_xy[1], z_mid],
                       [cw_xy[0], cw_xy[1], z_mid]])
    T_hw, T_cw = [], []
    for _, f in steps:
        m = pv.read(f)
        p = pv.PolyData(probes).sample(m)
        T_hw.append(float(p["T"][0])); T_cw.append(float(p["T"][1]))
    T_hw = np.array(T_hw); T_cw = np.array(T_cw)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(times_d, T_hw - 273.15, lw=2, color="tab:red", label="Hot Well")
    ax.plot(times_d, T_cw - 273.15, lw=2, color="tab:blue", label="Cold Well")
    ax.set_xlabel("Zeit [d]"); ax.set_ylabel("T [°C]")
    ax.set_title("Temperatur an HW und CW (Aquifer-Mitte)")
    ax.legend(); ax.grid(True, alpha=0.3); fig.tight_layout()
    fig.savefig(fig_dir / "1_well_temperature.png", dpi=130); plt.close(fig)

    # 2) Snapshots (xy-Slice durch Aquifer-Mitte) ----------------------
    idxs = np.linspace(0, len(steps) - 1, N_SNAPSHOTS).astype(int)
    fig, axes = plt.subplots(1, N_SNAPSHOTS, figsize=(4 * N_SNAPSHOTS, 4.5),
                              sharey=True)
    for ax, i in zip(axes, idxs):
        t, f = steps[i]
        m = pv.read(f)
        sl = m.slice(normal="z", origin=(0, 0, z_mid))
        pts = sl.points
        sc = ax.tricontourf(pts[:, 0], pts[:, 1], sl["T"] - 273.15,
                            levels=20, cmap="inferno")
        ax.scatter(*hw_xy, c="red",  s=40, marker="o", edgecolors="white", label="HW")
        ax.scatter(*cw_xy, c="blue", s=40, marker="o", edgecolors="white", label="CW")
        ax.set_xlabel("x [m]"); ax.set_title(f"t = {t/DAY:.0f} d"); ax.set_aspect("equal")
    axes[0].set_ylabel("y [m]")
    axes[-1].legend(loc="upper right")
    fig.colorbar(sc, ax=axes, shrink=0.7, label="T [°C]")
    fig.savefig(fig_dir / "2_field_snapshots.png", dpi=130); plt.close(fig)

    # 3) Energiebilanz -------------------------------------------------
    energy = []
    for _, f in steps:
        m = pv.read(f); dT = np.asarray(m["T"]) - T0
        energy.append(rho_cp * cartesian_integral(m, dT))
    energy = np.array(energy) / 1e9
    E_max = energy.max(); E_end = energy[-1]
    recovery = (E_max - E_end) / E_max if E_max > 0 else float("nan")

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(times_d, energy, lw=2)
    ax.axhline(E_max, ls="--", color="r", alpha=0.5, label=f"max. gespeichert: {E_max:.2f} GJ")
    ax.set_xlabel("Zeit [d]"); ax.set_ylabel("ΔE [GJ]")
    ax.set_title(f"Energiebilanz – Recovery ≈ {recovery*100:.0f} %")
    ax.legend(); ax.grid(True, alpha=0.3); fig.tight_layout()
    fig.savefig(fig_dir / "3_energy_balance.png", dpi=130); plt.close(fig)

    # 4) Plume-Reichweite um HW ---------------------------------------
    r_front = []
    for _, f in steps:
        m = pv.read(f); pts = m.points; dT = np.asarray(m["T"]) - T0
        in_aq = (pts[:, 2] >= z_aq_bot) & (pts[:, 2] <= z_aq_top)
        mask = in_aq & (dT > DELTA_T_THRESHOLD)
        if mask.any():
            d = np.hypot(pts[mask, 0] - hw_xy[0], pts[mask, 1] - hw_xy[1])
            r_front.append(float(d.max()))
        else:
            r_front.append(0.0)
    r_front = np.array(r_front)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(times_d, r_front, lw=2)
    ax.set_xlabel("Zeit [d]"); ax.set_ylabel(f"r(ΔT > {DELTA_T_THRESHOLD} K) um HW [m]")
    ax.set_title("Reichweite der warmen Fahne um den Hot Well")
    ax.grid(True, alpha=0.3); fig.tight_layout()
    fig.savefig(fig_dir / "4_plume_extent.png", dpi=130); plt.close(fig)

    print(f"Plots → {fig_dir}")
    print(f"  Recovery: {recovery*100:.1f} %  |  Emax: {E_max:.2f} GJ  |  max. Plume: {r_front.max():.1f} m")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
