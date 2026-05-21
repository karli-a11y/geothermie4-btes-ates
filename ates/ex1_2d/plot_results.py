#!/usr/bin/env python3
"""
Auswertung & Plots für die ATES-2D-Übung (radialsymmetrischer Single-Well).

Erzeugt in figures/:
  - 1_well_temperature.png       T(t) am Brunnenrand
  - 2_field_snapshots.png        T-Feld in (r,z) zu 4 Zeitpunkten
  - 3_energy_balance.png         Gespeicherte Wärmemenge & Recovery-Effizienz
  - 4_plume_extent.png           Reichweite der thermischen Front (ΔT > 1 K) im Aquifer

Anpassungen erfolgen im CONFIG des ates_radial_2d.py;
bei Bedarf DELTA_T_THRESHOLD / PROBE_R_OFFSET_M weiter unten ändern.
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

from ates_radial_2d import CONFIG

PROBE_R_OFFSET_M  = 0.5
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


def radial_integral(mesh, field):
    centers = mesh.cell_centers().points
    r = centers[:, 0]
    vols = mesh.compute_cell_sizes(length=False, area=True, volume=False)["Area"]
    if field.shape[0] == mesh.n_points:
        cell_field = np.empty(mesh.n_cells)
        for i in range(mesh.n_cells):
            pts = mesh.get_cell(i).point_ids
            cell_field[i] = field[pts].mean()
    else:
        cell_field = field
    return float(np.sum(cell_field * 2 * np.pi * r * vols))


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

    r_well = CONFIG["well"]["r_well_m"]
    z_aq_top = CONFIG["layers"]["caprock_bottom_thickness_m"] + CONFIG["layers"]["aquifer_thickness_m"]
    z_aq_bot = CONFIG["layers"]["caprock_bottom_thickness_m"]
    z_mid    = 0.5 * (z_aq_top + z_aq_bot)

    steps = read_pvd(pvd)
    times = np.array([t for t, _ in steps]); times_d = times / DAY

    # 1) T(t) ----------------------------------------------------------
    probe_point = np.array([[r_well + PROBE_R_OFFSET_M, z_mid, 0.0]])
    probe_T = []
    for _, f in steps:
        m = pv.read(f); probe_T.append(float(pv.PolyData(probe_point).sample(m)["T"][0]))
    probe_T = np.array(probe_T)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(times_d, probe_T - 273.15, lw=2)
    ax.set_xlabel("Zeit [d]"); ax.set_ylabel("T am Brunnen [°C]")
    ax.set_title(f"Temperatur am Brunnenrand (r={r_well + PROBE_R_OFFSET_M:.1f} m, Aquifer-Mitte)")
    ax.grid(True, alpha=0.3); fig.tight_layout()
    fig.savefig(fig_dir / "1_well_temperature.png", dpi=130); plt.close(fig)

    # 2) Snapshots — y-Achse als Tiefe unter Oberfläche -----------
    z_total = z_aq_top + CONFIG["layers"]["caprock_top_thickness_m"]
    idxs = np.linspace(0, len(steps) - 1, N_SNAPSHOTS).astype(int)
    fig, axes = plt.subplots(1, N_SNAPSHOTS, figsize=(4 * N_SNAPSHOTS, 5), sharey=True)
    for ax, i in zip(axes, idxs):
        t, f = steps[i]; m = pv.read(f); pts = m.points
        depth = z_total - pts[:, 1]
        sc = ax.tricontourf(pts[:, 0], depth, m["T"] - 273.15, levels=20, cmap="inferno")
        ax.axhline(z_total - z_aq_top, color="cyan", lw=0.6, ls="--")
        ax.axhline(z_total - z_aq_bot, color="cyan", lw=0.6, ls="--")
        ax.set_xlabel("r [m]"); ax.set_title(f"t = {t/DAY:.0f} d"); ax.set_aspect("equal")
        ax.invert_yaxis()
    axes[0].set_ylabel("Tiefe unter Oberfläche [m]")
    fig.colorbar(sc, ax=axes, shrink=0.7, label="T [°C]")
    fig.savefig(fig_dir / "2_field_snapshots.png", dpi=130); plt.close(fig)

    # 3) Energie -------------------------------------------------------
    energy = []
    for _, f in steps:
        m = pv.read(f); dT = np.asarray(m["T"]) - T0
        energy.append(rho_cp * radial_integral(m, dT))
    energy = np.array(energy) / 1e9
    E_max = energy.max(); E_end = energy[-1]
    # Theoretische Injektionsmenge: m_dot · c_p,f · ΔT · t_charge
    m_nom    = CONFIG["operation"]["mass_flow_rate_kg_s"]
    cp_f     = CONFIG["fluid"]["cp_J_kgK"]
    T_inj    = CONFIG["operation"]["T_hot_K"]
    dT_inj   = T_inj - T0
    t_charge = CONFIG["cycles"].get("charge_days", 91.25) * DAY
    E_inj_J  = m_nom * cp_f * dT_inj * t_charge
    recovery = (E_max - E_end) / E_max if E_max > 0 else float("nan")

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(times_d, energy, lw=2)
    ax.axhline(E_max, ls="--", color="r", alpha=0.5, label=f"max. gespeichert: {E_max:.2f} GJ")
    ax.set_xlabel("Zeit [d]"); ax.set_ylabel("ΔE [GJ]")
    ax.set_title(f"Energiebilanz – Recovery ≈ {recovery*100:.0f} %")
    ax.legend(); ax.grid(True, alpha=0.3); fig.tight_layout()
    fig.savefig(fig_dir / "3_energy_balance.png", dpi=130); plt.close(fig)

    # 4) Plume-Reichweite im Aquifer ----------------------------------
    r_front = []
    for _, f in steps:
        m = pv.read(f); pts = m.points; dT = np.asarray(m["T"]) - T0
        in_aq = (pts[:, 1] >= z_aq_bot) & (pts[:, 1] <= z_aq_top)
        mask = in_aq & (dT > DELTA_T_THRESHOLD)
        r_front.append(float(pts[mask, 0].max()) if mask.any() else 0.0)
    r_front = np.array(r_front)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(times_d, r_front, lw=2)
    ax.set_xlabel("Zeit [d]"); ax.set_ylabel(f"r(ΔT > {DELTA_T_THRESHOLD} K) im Aquifer [m]")
    ax.set_title("Reichweite der thermischen Front im Aquifer")
    ax.grid(True, alpha=0.3); fig.tight_layout()
    fig.savefig(fig_dir / "4_plume_extent.png", dpi=130); plt.close(fig)

    # --- 5) Massenstrom + Vorlauf­temperatur über Zeit -----------
    from ates_radial_2d import build_cycle_curves
    cv = build_cycle_curves(CONFIG)
    t_m, q_m = cv["cycle_mass"]
    t_T, q_T = cv["cycle_T"]
    m_nom = CONFIG["operation"]["mass_flow_rate_kg_s"]
    T_hot = CONFIG["operation"]["T_hot_K"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 5.5), sharex=True)
    ax1.plot(t_m / DAY, q_m * m_nom, lw=1.5, color="tab:orange")
    ax1.axhline(0, color="k", lw=0.6)
    ax1.set_ylabel("Massenstrom m_dot [kg/s]")
    ax1.set_title("Schaltprofil — Massenstrom und Vorlauf­temperatur")
    ax1.grid(True, alpha=0.3)
    ax2.plot(t_T / DAY, q_T * T_hot - 273.15, lw=1.5, color="tab:red")
    ax2.set_xlabel("Zeit [d]"); ax2.set_ylabel("Vorlauf-T [°C]")
    ax2.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "5_power_schedule.png", dpi=130); plt.close(fig)

    print(f"Plots → {fig_dir}")
    print(f"  Injektion pro Zyklus (m_dot·c_p·ΔT·t):  {E_inj_J/1e9:.2f} GJ")
    print(f"  max. integrierte Energie E_max:         {E_max:.2e} GJ")
    print(f"  Hinweis: E_max integriert über das ganze Domain inkl.")
    print(f"           durch die Dirichlet-T-BC implizit zugeführte Wärme")
    print(f"           — kann die reine Injektion deutlich übersteigen.")
    print(f"  Recovery η:                             {recovery*100:.1f} %")
    print(f"  max. Plume (ΔT > {DELTA_T_THRESHOLD} K):                {r_front.max():.1f} m")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
