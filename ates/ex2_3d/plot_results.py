#!/usr/bin/env python3
"""
Auswertung & Plots für die ATES-3D-Übung (Single-Well mit Grundwasserströmung).

Erzeugt in figures/:
  - 1_well_temperature.png       Brunnentemperatur (= Injektions-/dynamische
                                 Entnahmetemperatur) über den gesamten Zeitraum
  - 2_field_snapshots.png        T-Feld in der Aquifer-Mittelebene zu N_SNAPSHOTS
                                 Zeitpunkten (zeigt saisonales Atmen + GW-Drift)
  - 3_energy_balance.png         Gespeicherte Wärmemenge & Recovery-Effizienz
  - 4_plume_extent.png           Reichweite der thermischen Front um den Brunnen
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
N_SNAPSHOTS       = 8       # Feld-Snapshots über den Zeitraum (Raster N_ROWS×…)
N_ROWS            = 2       # Zeilen im Snapshot-Raster
VIEW_M            = 70.0    # halbe Ausschnittsbreite der Snapshots [m]
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
    z_aq_bot = CONFIG["layers"]["caprock_bottom_thickness_m"]
    z_aq_top = z_aq_bot + CONFIG["layers"]["aquifer_thickness_m"]
    z_mid    = 0.5 * (z_aq_bot + z_aq_top)

    steps = read_pvd(pvd)
    times = np.array([t for t, _ in steps]); times_d = times / DAY

    # --- Einzeldurchlauf: Brunnentemperatur, Energie, Fahnenreichweite ---
    #     (jede Datei nur EINMAL lesen)
    probe = pv.PolyData(np.array([[hw_xy[0], hw_xy[1], z_mid]]))
    T_hw    = np.empty(len(steps))
    energy  = np.empty(len(steps))
    r_front = np.empty(len(steps))
    for k, (_, f) in enumerate(steps):
        m = pv.read(f)
        T_hw[k] = float(probe.sample(m)["T"][0])
        dT = np.asarray(m["T"]) - T0
        energy[k] = rho_cp * cartesian_integral(m, dT) / 1e9
        pts = m.points
        in_aq = (pts[:, 2] >= z_aq_bot) & (pts[:, 2] <= z_aq_top)
        mask = in_aq & (dT > DELTA_T_THRESHOLD)
        r_front[k] = (float(np.hypot(pts[mask, 0] - hw_xy[0], pts[mask, 1] - hw_xy[1]).max())
                      if mask.any() else 0.0)

    # 1) Brunnentemperatur über die Zeit ------------------------------

    T_amb_c = T0 - 273.15; T_inj_c = CONFIG["operation"]["T_hot_K"] - 273.15
    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(times_d, T_hw - 273.15, lw=1.8, color="#b22222", label="Brunnentemperatur")
    ax.axhline(T_inj_c, color="k", lw=0.8, ls="--", label=f"T_inj = {T_inj_c:.0f} °C")
    ax.axhline(T_amb_c, color="k", lw=0.8, ls=":",  label=f"T_amb = {T_amb_c:.0f} °C")
    ax.set_xlabel("Zeit [d]"); ax.set_ylabel("Brunnentemperatur [°C]")
    ax.set_title("ATES 3D — Brunnentemperatur: Injektion bei Beladung, dynamische Entnahme bei Förderung")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3); fig.tight_layout()
    fig.savefig(fig_dir / "1_well_temperature.png", dpi=130); plt.close(fig)

    # 2) Snapshots (xy-Slice durch Aquifer-Mitte) — N_SNAPSHOTS Zeitpunkte ----
    from matplotlib.tri import Triangulation
    T_amb = T0 - 273.15; T_inj = CONFIG["operation"]["T_hot_K"] - 273.15
    gw_on = CONFIG.get("regional_gw", {}).get("enable", False)
    idxs = np.unique(np.linspace(0, len(steps) - 1, N_SNAPSHOTS).astype(int))
    ncol = int(np.ceil(len(idxs) / N_ROWS))
    fig, axes = plt.subplots(N_ROWS, ncol, figsize=(3.6 * ncol, 3.6 * N_ROWS),
                             sharex=True, sharey=True, squeeze=False)
    axes = axes.ravel()
    levels = np.linspace(T_amb, T_inj, 26); sc = None
    for ax, i in zip(axes, idxs):
        t, f = steps[i]
        sl = pv.read(f).slice(normal="z", origin=(0, 0, z_mid))
        pts = sl.points; tri = Triangulation(pts[:, 0], pts[:, 1])
        sc = ax.tricontourf(tri, sl["T"] - 273.15, levels=levels,
                            cmap="coolwarm", extend="both")
        ax.tricontour(tri, sl["T"] - 273.15,
                      levels=[T_amb + 5, 0.5 * (T_amb + T_inj)], colors="k", linewidths=0.3)
        ax.plot(hw_xy[0], hw_xy[1], "ko", ms=4)
        ax.set_title(f"t = {t/DAY:.0f} d", fontsize=9)
        ax.set_xlim(-VIEW_M, VIEW_M); ax.set_ylim(-VIEW_M, VIEW_M); ax.set_aspect("equal")
    for ax in axes[len(idxs):]:
        ax.axis("off")
    if gw_on:
        axes[0].annotate("GW →", (0.66, 0.88), xycoords="axes fraction",
                         fontsize=10, weight="bold")
    fig.suptitle("ATES 3D — T-Feld in der Aquifer-Mitte: saisonales Atmen + GW-Drift",
                 fontsize=12)
    if sc is not None:
        fig.colorbar(sc, ax=list(axes), shrink=0.85, label="T [°C]")
    fig.savefig(fig_dir / "2_field_snapshots.png", dpi=130, bbox_inches="tight"); plt.close(fig)

    # 3) Energiebilanz (energy bereits im Einzeldurchlauf berechnet) ---
    E_max = energy.max(); E_end = energy[-1]
    recovery = (E_max - E_end) / E_max if E_max > 0 else float("nan")

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(times_d, energy, lw=2)
    ax.axhline(E_max, ls="--", color="r", alpha=0.5, label=f"max. gespeichert: {E_max:.2f} GJ")
    ax.set_xlabel("Zeit [d]"); ax.set_ylabel("ΔE [GJ]")
    ax.set_title(f"Energiebilanz – Recovery ≈ {recovery*100:.0f} %")
    ax.legend(); ax.grid(True, alpha=0.3); fig.tight_layout()
    fig.savefig(fig_dir / "3_energy_balance.png", dpi=130); plt.close(fig)

    # 4) Fahnenreichweite (r_front bereits im Einzeldurchlauf berechnet) --
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
