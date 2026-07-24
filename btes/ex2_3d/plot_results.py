#!/usr/bin/env python3
"""
Auswertung & Plots für die BTES-3D-Übung (Sondenfeld, konduktiver Wärmeeintrag).

Erzeugt in figures/:
  - 1_feasibility.png            MACHBARKEIT: Sondentemperatur (heißeste/kälteste
                                 Sonde) vs. Betriebsgrenzen + spezifische Rate W/m
  - 2_field_snapshots.png        T-Feld in der Feld-Mittelebene zu N_SNAPSHOTS
                                 Zeitpunkten (saisonales Atmen + Feld-Interaktion)
  - 3_energy_balance.png         Gespeicherte Wärmemenge & Recovery-Effizienz
  - 4_plume_extent.png           Radius der thermischen Front (ΔT > 1 K)

Machbarkeit (deine Frage): Bei einem DEFINIERTEN Wärmeeintrag (W bzw. W/m) ist
das Feld realistisch, wenn die Sondenwand-Temperatur in den Betriebsgrenzen
bleibt (Beladung ≤ T_charge_max, Förderung ≥ T_discharge_min). Läuft sie darüber
hinaus → Feld zu klein/eng oder Leistung zu hoch.
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

DELTA_T_THRESHOLD = 1.0
N_SNAPSHOTS       = 8       # Feld-Snapshots über den Zeitraum
N_ROWS            = 2       # Raster N_ROWS × …
VIEW_M            = 30.0    # halbe Ausschnittsbreite der Snapshots [m]
DAY = 86400.0


def load_config():
    """Erkennt automatisch, welche BTES-3D-Variante gelaufen ist, und lädt CONFIG."""
    here = Path(__file__).parent
    for prefix, module in (("btes_3d_bhe", "btes_3d_bhe"), ("btes_3d", "btes_3d")):
        if (here / "out" / f"{prefix}.pvd").exists():
            return __import__(module).CONFIG
    return __import__("btes_3d").CONFIG


def temperature_field(mesh) -> str:
    for name in ("temperature_soil", "T"):
        if name in mesh.point_data:
            return name
    raise KeyError("Kein Temperaturfeld (temperature_soil/T) im Mesh gefunden.")


def read_pvd(pvd_path: Path):
    tree = ET.parse(pvd_path); items = []
    base = pvd_path.parent
    for ds in tree.iter("DataSet"):
        items.append((float(ds.attrib["timestep"]), base / ds.attrib["file"]))
    items.sort(); return items


def effective_rho_cp(cfg):
    layers = cfg["layers"]
    soil = layers[len(layers) // 2]
    phi = soil["porosity"]
    return (phi * cfg["fluid"]["rho_ref_kg_m3"] * cfg["fluid"]["cp_J_kgK"]
            + (1 - phi) * soil["rho_s_kg_m3"] * soil["cp_s_J_kgK"])


def borehole_positions(cfg):
    fld = cfg["field"]
    if fld.get("positions"):
        return np.array(fld["positions"])
    nx, ny, sp = fld["n_x"], fld["n_y"], fld["spacing_m"]
    xs = (np.arange(nx) - (nx - 1) / 2) * sp
    ys = (np.arange(ny) - (ny - 1) / 2) * sp
    return np.array([(x, y) for y in ys for x in xs])


def cartesian_integral(mesh, field):
    vols = mesh.compute_cell_sizes(length=False, area=False, volume=True)["Volume"]
    if field.shape[0] == mesh.n_points:
        # vektorisiert: Punkt- -> Zelldaten (statt Python-Schleife über Zellen)
        tmp = mesh.copy(deep=False)
        tmp.point_data["_intg"] = np.asarray(field)
        cell_field = np.asarray(tmp.point_data_to_cell_data().cell_data["_intg"])
    else:
        cell_field = field
    return float(np.sum(cell_field * vols))


def main() -> int:
    CONFIG = load_config()
    out_dir = Path(CONFIG["output"]["out_dir"]); prefix = CONFIG["output"]["prefix"]
    pvd = out_dir / f"{prefix}.pvd"
    if not pvd.exists():
        print(f"FEHLER: {pvd} nicht gefunden. Bitte zuerst die Simulation ausführen.")
        return 1
    fig_dir = Path(__file__).parent / "figures"; fig_dir.mkdir(exist_ok=True)

    T0 = CONFIG["initial"]["T_K"]; T0_c = T0 - 273.15
    rho_cp = effective_rho_cp(CONFIG)
    positions = borehole_positions(CONFIG)

    # Geometrie: z=0 unten → Oberfläche bei z_total; Sondentiefe von oben gemessen
    z_total = sum(L["thickness_m"] for L in CONFIG["layers"])
    z_field_top = z_total - CONFIG["borehole"]["depth_top_m"]
    z_field_bot = z_total - CONFIG["borehole"]["depth_bottom_m"]
    z_mid = 0.5 * (z_field_top + z_field_bot)
    bh_len = CONFIG["borehole"]["depth_bottom_m"] - CONFIG["borehole"]["depth_top_m"]

    # Betriebsgrenzen + spezifische Rate (robust ggü. Config-Varianten)
    rl = CONFIG.get("realism", {})
    T_hi = rl.get("T_charge_max_C", 60.0); T_lo = rl.get("T_discharge_min_C", 4.0)
    P_ref = CONFIG.get("operation", {}).get("power_per_borehole_W", 0.0)
    mp = CONFIG.get("cycles", {}).get("monthly_power_W")
    P_peak = max(abs(P_ref), max((abs(p) for p in mp), default=0.0) if mp else 0.0)
    spec_rate = P_peak / bh_len   # W/m

    steps = read_pvd(pvd)
    times_d = np.array([t for t, _ in steps]) / DAY
    TEMP = temperature_field(pv.read(steps[0][1]))

    # --- Einzeldurchlauf: Sondentemperaturen, Energie, Fahnenradius ---
    probe = pv.PolyData(np.column_stack([positions, np.full(len(positions), z_mid)]))
    T_bh = np.empty((len(steps), len(positions)))    # je Sonde
    energy = np.empty(len(steps)); r_front = np.empty(len(steps))
    cx, cy = positions.mean(axis=0)
    for k, (_, f) in enumerate(steps):
        m = pv.read(f)
        T_bh[k] = np.asarray(probe.sample(m)[TEMP])
        dT = np.asarray(m[TEMP]) - T0
        energy[k] = rho_cp * cartesian_integral(m, dT) / 1e9
        pts = m.points
        in_zone = (pts[:, 2] >= z_field_bot) & (pts[:, 2] <= z_field_top)
        mask = in_zone & (dT > DELTA_T_THRESHOLD)
        r_front[k] = (float(np.hypot(pts[mask, 0] - cx, pts[mask, 1] - cy).max())
                      if mask.any() else 0.0)
    T_bh_c = T_bh - 273.15
    T_max = T_bh_c.max(axis=1)   # heißeste Sonde (Feldmitte bei Beladung)
    T_min = T_bh_c.min(axis=1)   # kälteste Sonde (bei Förderung)

    # 1) MACHBARKEIT: Sondentemperatur vs. Betriebsgrenzen ------------
    over_hi = T_max.max() > T_hi; under_lo = T_min.min() < T_lo
    feasible = not (over_hi or under_lo)
    fig, ax = plt.subplots(figsize=(12, 4.6))
    ax.fill_between(times_d, T_hi, max(T_hi + 5, T_max.max() + 3), color="red", alpha=0.06)
    ax.fill_between(times_d, min(T_lo - 5, T_min.min() - 3), T_lo, color="blue", alpha=0.06)
    ax.plot(times_d, T_max, color="#c0392b", lw=1.8, label="heißeste Sonde (Feldmitte)")
    ax.plot(times_d, T_min, color="#2471a3", lw=1.8, label="kälteste Sonde")
    ax.axhline(T_hi, color="red",  lw=1.0, ls="--", label=f"Grenze Beladung {T_hi:.0f} °C")
    ax.axhline(T_lo, color="blue", lw=1.0, ls="--", label=f"Grenze Förderung {T_lo:.0f} °C")
    ax.axhline(T0_c, color="k", lw=0.6, ls=":")
    ax.set_xlabel("Zeit [d]"); ax.set_ylabel("Sondentemperatur [°C]")
    verdict = "MACHBAR ✓" if feasible else "NICHT MACHBAR ✗ (Grenze überschritten)"
    spec_ok = rl.get("spec_rate_min_W_m", 20) <= spec_rate <= rl.get("spec_rate_max_W_m", 70)
    ax.set_title(f"BTES-Machbarkeit — spez. Rate {spec_rate:.0f} W/m "
                 f"({'im' if spec_ok else 'außerhalb'} Literaturbereich {rl.get('spec_rate_min_W_m',20):.0f}–{rl.get('spec_rate_max_W_m',70):.0f}) — {verdict}")
    ax.legend(loc="upper right", fontsize=8, ncol=2); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(fig_dir / "1_feasibility.png", dpi=130); plt.close(fig)

    # 2) Snapshots (xy-Slice in Feld-Mitteltiefe) — N_SNAPSHOTS -------
    from matplotlib.tri import Triangulation
    idxs = np.unique(np.linspace(0, len(steps) - 1, N_SNAPSHOTS).astype(int))
    ncol = int(np.ceil(len(idxs) / N_ROWS))
    fig, axes = plt.subplots(N_ROWS, ncol, figsize=(3.5 * ncol, 3.5 * N_ROWS),
                             sharex=True, sharey=True, squeeze=False)
    axes = axes.ravel()
    # symmetrische Skala um T0 (coolwarm): blau = Förderung/kalt, rot = Beladung/warm
    dT_amp = max(T_max.max() - T0_c, T0_c - T_min.min(), 2.0)
    levels = np.linspace(T0_c - dT_amp, T0_c + dT_amp, 24); sc = None
    for ax, i in zip(axes, idxs):
        t, f = steps[i]
        sl = pv.read(f).slice(normal="z", origin=(0, 0, z_mid))
        pts = sl.points; tri = Triangulation(pts[:, 0], pts[:, 1])
        sc = ax.tricontourf(tri, sl[TEMP] - 273.15, levels=levels, cmap="coolwarm", extend="both")
        ax.plot(positions[:, 0], positions[:, 1], "k+", ms=6, mew=1.2)
        ax.set_title(f"t = {t/DAY:.0f} d", fontsize=9)
        ax.set_xlim(-VIEW_M, VIEW_M); ax.set_ylim(-VIEW_M, VIEW_M); ax.set_aspect("equal")
    for ax in axes[len(idxs):]:
        ax.axis("off")
    fig.suptitle("BTES 3D — T-Feld in Sonden-Mitteltiefe: saisonales Atmen + Feld-Interaktion", fontsize=12)
    if sc is not None:
        fig.colorbar(sc, ax=list(axes), shrink=0.85, label="T [°C]")
    fig.savefig(fig_dir / "2_field_snapshots.png", dpi=130, bbox_inches="tight"); plt.close(fig)

    # 3) Gespeicherte Wärme (Beladung +, Förderung −; schwingt um 0) ---
    #    (Ein klassischer Recovery-%-Wert ist bei BTES mit ΣP=0 nicht sinnvoll,
    #     da die Energie um 0 pendelt statt monoton entladen zu werden.)
    fig, ax = plt.subplots(figsize=(11, 4))
    ax.axhline(0.0, color="k", lw=0.6, ls=":")
    ax.fill_between(times_d, 0, energy, where=(energy >= 0), color="red",  alpha=0.15)
    ax.fill_between(times_d, 0, energy, where=(energy < 0),  color="blue", alpha=0.15)
    ax.plot(times_d, energy, "k-", lw=1.6)
    ax.set_xlabel("Zeit [d]"); ax.set_ylabel("gespeicherte Wärme über T0 [GJ]")
    ax.set_title(f"BTES — gespeicherte Wärme im Untergrund "
                 f"(Beladung max +{energy.max():.1f} GJ, Förderung min {energy.min():.1f} GJ)")
    ax.grid(alpha=0.3); fig.tight_layout()
    fig.savefig(fig_dir / "3_energy_balance.png", dpi=130); plt.close(fig)

    # 4) Fahnenreichweite ---------------------------------------------
    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(times_d, r_front, lw=1.8)
    ax.set_xlabel("Zeit [d]"); ax.set_ylabel(f"r(ΔT > {DELTA_T_THRESHOLD} K) [m]")
    ax.set_title("BTES — Reichweite der thermischen Front (horizontaler Radius vom Feldzentrum)")
    ax.grid(alpha=0.3); fig.tight_layout()
    fig.savefig(fig_dir / "4_plume_extent.png", dpi=130); plt.close(fig)

    print(f"Plots → {fig_dir}")
    print(f"  Spez. Rate: {spec_rate:.1f} W/m  |  T_Sonde: {T_min.min():.1f}…{T_max.max():.1f} °C"
          f"  |  Grenzen {T_lo:.0f}/{T_hi:.0f} °C -> {'MACHBAR' if feasible else 'NICHT MACHBAR'}")
    print(f"  gespeicherte Wärme: {energy.min():.1f}…{energy.max():.1f} GJ  |  max. Fahne: {r_front.max():.1f} m")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
