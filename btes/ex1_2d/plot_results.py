#!/usr/bin/env python3
"""
Auswertung & Plots für die BTES-2D-Übung.

Liest die OGS-Ausgaben (out/*.pvd, *.vtu) und erzeugt in figures/:
  - 1_well_temperature.png       T(t) am Sondenrand
  - 2_field_snapshots.png        T-Feld in (r,z) zu 4 Zeitpunkten
  - 3_energy_balance.png         Gespeicherte Wärmemenge & Recovery-Effizienz
  - 4_plume_extent.png           Reichweite der thermischen Front (ΔT > 1 K)

Hier nichts ändern – ggf. nur PROBE_R_M / DELTA_T_THRESHOLD
unten anpassen, wenn nötig.
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

from btes_radial_2d import CONFIG

# ---- Stellschrauben (selten nötig) ----------------------------------
PROBE_R_M          = 0.5      # Auswertepunkt: r = Sondenrand
DELTA_T_THRESHOLD  = 1.0      # Schwelle für „thermische Front"
N_SNAPSHOTS        = 4
# ----------------------------------------------------------------------

DAY = 86400.0


def read_pvd(pvd_path: Path) -> list[tuple[float, Path]]:
    tree = ET.parse(pvd_path)
    items = []
    base = pvd_path.parent
    for ds in tree.iter("DataSet"):
        t = float(ds.attrib["timestep"])
        f = base / ds.attrib["file"]
        items.append((t, f))
    items.sort()
    return items


def effective_rho_cp(cfg: dict) -> float:
    soil = cfg["materials"]["soil"]
    phi  = soil["porosity"]
    rho_s = soil["rho_s_kg_m3"]
    cp_s  = soil["cp_s_J_kgK"]
    rho_f = cfg["fluid"]["rho_ref_kg_m3"]
    cp_f  = cfg["fluid"]["cp_J_kgK"]
    return phi * rho_f * cp_f + (1.0 - phi) * rho_s * cp_s


def radial_integral(mesh: pv.UnstructuredGrid, field: np.ndarray) -> float:
    """∫ field · 2π r dV über das radial-2D-Mesh (Achsendrehung um r=0)."""
    centers = mesh.cell_centers().points
    r = centers[:, 0]
    vols = mesh.compute_cell_sizes(length=False, area=True, volume=False)["Area"]
    # Punktdaten → Zellmittel
    if field.shape[0] == mesh.n_points:
        cdata = mesh.point_data_to_cell_data().point_data
        # Fallback: einfacher Mittelwert über Zellpunkte
        cell_field = np.empty(mesh.n_cells)
        for i in range(mesh.n_cells):
            pts = mesh.get_cell(i).point_ids
            cell_field[i] = field[pts].mean()
    else:
        cell_field = field
    return float(np.sum(cell_field * 2.0 * np.pi * r * vols))


def main() -> int:
    out_dir = Path(CONFIG["output"]["out_dir"])
    prefix  = CONFIG["output"]["prefix"]
    pvd     = out_dir / f"{prefix}.pvd"
    if not pvd.exists():
        print(f"FEHLER: {pvd} nicht gefunden. Bitte zuerst die Simulation ausführen:")
        print(f"        python {Path(__file__).parent.name}/btes_radial_2d.py")
        return 1

    fig_dir = Path(__file__).parent / "figures"
    fig_dir.mkdir(exist_ok=True)

    T0 = CONFIG["initial"]["T_K"]
    rho_cp = effective_rho_cp(CONFIG)

    steps = read_pvd(pvd)
    times = np.array([t for t, _ in steps])
    times_d = times / DAY

    # --- 1) T(t) am Sondenrand --------------------------------------
    # Gesamttiefe = Summe der Schicht­dicken. Sondenmitte in z (internal):
    #   z_top liegt an der Oberfläche, z = 0 am Modellboden.
    z_total = sum(L["thickness_m"] for L in CONFIG["layers"])
    bh = CONFIG["borehole"]
    z_bh_mid = z_total - 0.5 * (bh["depth_top_m"] + bh["depth_bottom_m"])
    probe_T = []
    probe_point = np.array([[PROBE_R_M, z_bh_mid, 0.0]])
    for _, f in steps:
        m = pv.read(f)
        probed = pv.PolyData(probe_point).sample(m)
        probe_T.append(float(probed["T"][0]))
    probe_T = np.array(probe_T)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(times_d, probe_T - 273.15, lw=2)
    ax.set_xlabel("Zeit [d]")
    ax.set_ylabel(f"T am Sondenrand (r={PROBE_R_M} m) [°C]")
    ax.set_title("Temperatur am Sondenrand")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "1_well_temperature.png", dpi=130)
    plt.close(fig)

    # --- 2) Feld-Snapshots ------------------------------------------
    # y-Achse als Tiefe unter Oberfläche dargestellt (depth = z_total - z).
    idxs = np.linspace(0, len(steps) - 1, N_SNAPSHOTS).astype(int)
    fig, axes = plt.subplots(1, N_SNAPSHOTS, figsize=(4 * N_SNAPSHOTS, 5),
                              sharey=True)
    for ax, i in zip(axes, idxs):
        t, f = steps[i]
        m = pv.read(f)
        pts = m.points
        T = m["T"]
        depth = z_total - pts[:, 1]
        sc = ax.tricontourf(pts[:, 0], depth, T - 273.15,
                            levels=20, cmap="inferno")
        # Schicht­grenzen einzeichnen
        z_acc = 0.0
        for L in CONFIG["layers"][:-1]:
            z_acc += L["thickness_m"]
            ax.axhline(z_acc, color="gray", lw=0.5, ls="--")
        ax.set_xlabel("r [m]")
        ax.set_title(f"t = {t/DAY:.0f} d")
        ax.set_aspect("equal")
        ax.invert_yaxis()
    axes[0].set_ylabel("Tiefe unter Oberfläche [m]")
    fig.colorbar(sc, ax=axes, shrink=0.7, label="T [°C]")
    fig.savefig(fig_dir / "2_field_snapshots.png", dpi=130)
    plt.close(fig)

    # --- 3) Energiebilanz & Recovery --------------------------------
    energy = []
    for _, f in steps:
        m = pv.read(f)
        dT = np.asarray(m["T"]) - T0
        E = rho_cp * radial_integral(m, dT)
        energy.append(E)
    energy = np.array(energy) / 1e9  # GJ

    E_charged = energy.max()
    E_end     = energy[-1]
    # Theoretische Injektionsmenge pro Lade­phase (Q_h * V * t_charge)
    P_ref     = CONFIG["operation"]["power_per_borehole_W"]
    t_charge  = CONFIG["cycles"].get("charge_days", 91.25) * DAY
    E_injected_per_cycle_J = P_ref * t_charge
    # Recovery-Effizienz: Anteil der gespeicherten Energie, der wieder
    # entzogen wird. Werte > 100 % treten auf bei numerischer Drift im
    # symmetrischen Q_h-Zyklus (Bilanz endet knapp unter null) — das
    # ist ein Modell-Artefakt und kein physikalisches Ergebnis.
    # Für realistische Verluste asymmetrische Lade-/Förder­dauer setzen
    # (siehe Aufgaben in der README).
    recovery = (E_charged - E_end) / E_charged if E_charged > 0 else float("nan")

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(times_d, energy, lw=2)
    ax.axhline(E_charged, ls="--", color="r", alpha=0.5,
               label=f"max. gespeichert: {E_charged:.2f} GJ")
    ax.set_xlabel("Zeit [d]")
    ax.set_ylabel("Gespeicherte Wärmemenge ΔE [GJ]")
    ax.set_title(f"Energiebilanz – Recovery ≈ {recovery*100:.0f} %")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "3_energy_balance.png", dpi=130)
    plt.close(fig)

    # --- 4) Plume-Reichweite ----------------------------------------
    r_front = []
    for _, f in steps:
        m = pv.read(f)
        pts = m.points
        dT = np.asarray(m["T"]) - T0
        mask = dT > DELTA_T_THRESHOLD
        r_front.append(float(pts[mask, 0].max()) if mask.any() else 0.0)
    r_front = np.array(r_front)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(times_d, r_front, lw=2)
    ax.set_xlabel("Zeit [d]")
    ax.set_ylabel(f"r(ΔT > {DELTA_T_THRESHOLD} K) [m]")
    ax.set_title("Reichweite der thermischen Front")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "4_plume_extent.png", dpi=130)
    plt.close(fig)

    # --- 5) Leistungs-/Zyklus-Schedule (inkl. Rampe) ----------------
    from btes_radial_2d import build_cycle_curves
    cv = build_cycle_curves(CONFIG)
    t_c, q_c = cv["cycle_q"]
    P_ref = CONFIG["operation"]["power_per_borehole_W"]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(t_c / DAY, q_c * P_ref / 1000.0, lw=1.5, color="tab:orange")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("Zeit [d]")
    ax.set_ylabel("Sondenleistung [kW]")
    ax.set_title("Schaltprofil P(t) — Lade (+) / Förder (−) inkl. Rampen")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "5_power_schedule.png", dpi=130)
    plt.close(fig)

    print(f"Plots geschrieben nach: {fig_dir}")
    print(f"  Injektion pro Zyklus (P·t_charge):  {E_injected_per_cycle_J/1e9:.2f} GJ")
    print(f"  max. gespeicherte Energie E_max:    {E_charged:.2f} GJ")
    print(f"  Recovery-Effizienz η:               {recovery*100:.1f} %")
    if recovery > 1.001:
        print(f"  ! η > 100 % deutet auf symmetrischen Zyklus + numerische Drift")
        print(f"    (E_end knapp negativ). Asym­metrische Phasen­dauer für")
        print(f"    realistische Verlust-Quantifizierung verwenden.")
    print(f"  max. Reichweite (ΔT > {DELTA_T_THRESHOLD} K):  {r_front.max():.1f} m")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
