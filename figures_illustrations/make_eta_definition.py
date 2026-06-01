"""Regenerate eta_definition.png — schematic recovery-efficiency balance.

The original illustration script was lost; this reproduces the schematic used
in HANDBUCH.md (Section 9.3) with two corrections:
  1. The stored-heat curve DeltaE ends at a positive residual E_end > 0
     (a residual always remains: conduction into cap/base layers, groundwater
     advection, incomplete recovery), so eta < 100 %. The residual is drawn as
     a variable/uncertain band rather than a single value.
  2. The E_end annotation/arrow is placed in the free white area to the right
     of the phase bands so its label no longer collides with the "Storage"
     phase label.

Run:  python make_eta_definition.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch  # noqa: F401 (style ref)

OUT = Path(__file__).with_name("eta_definition.png")

NAVY = "#1f3b6e"
ORANGE = "#d2691e"

# Phase boundaries [days]
t_charge_end, t_store1_end, t_disch_end, t_store2_end = 90, 182, 273, 365

# Residual at cycle end (variable / plant-dependent), drawn as a band
E_end = 0.13
E_end_lo, E_end_hi = 0.05, 0.22

fig, ax = plt.subplots(figsize=(10, 5))

# Phase background bands
bands = [
    (0, t_charge_end, "#fbe3d6", "Lade"),
    (t_charge_end, t_store1_end, "#eeeeee", "Storage"),
    (t_store1_end, t_disch_end, "#dde7f2", "Förder"),
    (t_disch_end, t_store2_end, "#eeeeee", "Storage"),
]
for x0, x1, colour, label in bands:
    ax.axvspan(x0, x1, color=colour, alpha=0.9, lw=0)
    ax.text((x0 + x1) / 2, -0.13, label, ha="center", va="center",
            fontsize=11, color="#555555")

# Reference lines
ax.axhline(1.0, ls="--", lw=1.2, color=ORANGE, alpha=0.8)
ax.axhline(0.0, ls="--", lw=1.0, color="#888888")

# Stored-heat curve DeltaE(t): rise -> slight storage loss -> discharge -> residual
tx = [0, t_charge_end, t_store1_end, t_disch_end, t_store2_end]
ey = [0.0, 1.0, 0.90, 0.18, E_end]
ax.plot(tx, ey, "-", color=NAVY, lw=3.2, solid_capstyle="round", zorder=5)

# E_max annotation (kept near the peak)
ax.annotate("E_max", xy=(t_charge_end, 1.0), xytext=(t_charge_end + 18, 1.10),
            color=ORANGE, fontsize=13, va="center",
            arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.5))

# Endpoint marker
ax.plot([t_store2_end], [E_end], "o", color=NAVY, ms=7, zorder=6)

# Variable residual band at the cycle end (free white area to the right)
x_band0, x_band1 = t_store2_end, 452
ax.fill_between([x_band0, x_band1], E_end_lo, E_end_hi,
                color=NAVY, alpha=0.12, lw=0)
ax.plot([x_band0, x_band1], [E_end, E_end], ls=":", lw=1.3, color=NAVY, alpha=0.7)

# E_end annotation in the clear white region (no overlap with phase labels)
ax.annotate(r"$E_\mathrm{end}>0$" "\n(Rest, anlagenabhängig)",
            xy=(t_store2_end, E_end), xytext=(395, 0.46),
            color=NAVY, fontsize=11.5, ha="center", va="center",
            arrowprops=dict(arrowstyle="->", color=NAVY, lw=1.5))

# Efficiency formula box (right, mid-height — clear of curve and labels)
ax.text(452, 0.78, r"$\eta = (E_\mathrm{max} - E_\mathrm{end})\,/\,E_\mathrm{max}$",
        fontsize=14, color=NAVY, ha="right", va="center",
        bbox=dict(boxstyle="round,pad=0.5", fc="white", ec=NAVY, lw=1.5))

ax.set_xlim(0, 455)
ax.set_ylim(-0.22, 1.22)
ax.set_xlabel("Zeit [d]", fontsize=12)
ax.set_ylabel("Gespeicherte Wärmemenge ΔE  (normiert)", fontsize=12)
ax.set_title("Recovery-Effizienz η — schematische Bilanz über einen Zyklus",
             fontsize=13)
ax.grid(True, axis="y", ls=":", lw=0.6, color="#cccccc", alpha=0.6)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)

fig.tight_layout()
fig.savefig(OUT, dpi=130)
print(f"wrote {OUT}")
