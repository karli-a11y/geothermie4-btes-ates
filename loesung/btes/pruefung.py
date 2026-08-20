#!/usr/bin/env python3
"""Selbstprüfung des Modells.

Rechnet eine Reihe kleiner Fälle und vergleicht sie gegen Größen, die aus der
Physik oder der Geometrie exakt folgen. Läuft in wenigen Minuten.

    python pruefung.py

Jede Prüfung meldet OK oder FEHLER samt gemessener Abweichung. Am Ende steht
eine Zusammenfassung; der Rückgabewert ist 0 nur, wenn alles hält.
"""
from __future__ import annotations

import copy
import glob
import re
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import btes_loesung as B

ERGEBNIS: list[tuple[str, bool, str]] = []


def pruefe(name, bedingung, detail=""):
    ERGEBNIS.append((name, bool(bedingung), detail))
    print(f"  [{'OK  ' if bedingung else 'FEHL'}] {name}"
          + (f"   {detail}" if detail else ""))
    return bool(bedingung)


def klein(nx=4, ny=4, jahre=1, symmetrie="viertel", zelle=False):
    """Winziger Testfall — schnell genug für eine Selbstprüfung."""
    c = copy.deepcopy(B.CONFIG)
    c["ablauf"] = {k: False for k in c["ablauf"]}
    c["field"].update(n_x_full=nx, n_y_full=ny, abstand_m=6.0,
                      symmetrie=symmetrie, einheitszelle=zelle)
    c["zeit"].update(jahre=jahre, schritte_je_monat=2, ausgabe_je_n_schritte=6)
    c["netz"].update(automatisch=True, ziel_knoten=400_000,
                     elemente_je_box=2.0, dz_je_sonde=12)
    c["domain"].update(automatisch=False, size_x_m=40.0, size_y_m=40.0)
    return c


def lauf(c, d: Path):
    d.mkdir(parents=True, exist_ok=True)
    c["ausgabe"] = dict(c["ausgabe"]); c["ausgabe"]["ordner"] = str(d)
    B.netz_bauen(c, d)
    nf = B.netz_wandeln(c, d / f"{c['ausgabe']['prefix']}.msh", d)
    rc = B.ogs_starten(B.prj_bauen(c, d, nf, B.lastkurve(c)))
    if rc != 0:
        raise RuntimeError(f"OGS endete mit {rc}")
    return d


def waermeinhalt(c, d: Path):
    """Wärmeinhalt am letzten Ausgabezeitpunkt und Volumina."""
    import pyvista as pv
    pref = c["ausgabe"]["prefix"]
    f = sorted(glob.glob(str(d / f"{pref}_ts_*.vtu")),
               key=lambda x: float(re.search(r"_t_([\d.]+)\.vtu$", x).group(1)))[-1]
    dom = pv.read(str(d / f"{pref}_domain.vtu"))
    mid = np.asarray(dom.cell_data["MaterialIDs"])
    vol = np.abs(np.asarray(dom.compute_cell_sizes(
        length=False, area=False, volume=True).cell_data["Volume"]))
    lb = list(reversed(c["layers"]))
    rc = np.zeros(dom.n_cells)
    for i, L in enumerate(lb):
        rc[mid == i] = B.effektive_stoffwerte(L, c["fluid"])["rho_c"]
    rc[mid >= len(lb)] = B.effektive_stoffwerte(
        lb[len(lb) // 2], c["fluid"])["rho_c"]
    T = np.asarray(pv.read(f).point_data_to_cell_data().cell_data["T"])
    T0 = c["initial"]["T_C"] + 273.15
    return (float((rc * (T - T0) * vol).sum()) / 1e9, float(vol.sum()),
            float(vol[mid >= len(lb)].sum()))


# ======================================================================
print("=" * 70)
print("1  Aufbau und Lastrechnung")
print("=" * 70)

c = copy.deepcopy(B.CONFIG)
pruefe("CONFIG hat die erwarteten Blöcke",
       set(c) >= {"ablauf", "last", "field", "zeit", "borehole", "domain",
                  "netz", "layers", "fluid", "initial", "ausgabe", "loeser"})

lp0 = B.lastprofil(c)
pruefe("Reduktionsgrad 0 reproduziert die Aufgabenstellung",
       abs(lp0["bedarf_MWh"] - c["last"]["waermebedarf_MWh_a"]) < 1,
       f"{lp0['bedarf_MWh']:,.0f} MWh/a")

r_null = lp0["reduktion_bilanz_null"]
c2 = copy.deepcopy(c); c2["last"]["reduktionsgrad_prozent"] = r_null
lpb = B.lastprofil(c2)
pruefe("gemeldeter Reduktionsgrad schließt die Bilanz",
       abs(lpb["bilanz_MWh"]) < 1.0,
       f"{r_null:.2f} % -> Bilanz {lpb['bilanz_MWh']:+.2f} MWh/a")
pruefe("bei geschlossener Bilanz ist Entnahme = Beladung",
       abs(lpb["entnahme_MWh"] - lpb["beladung_MWh"]) < 1.0,
       f"{lpb['entnahme_MWh']:,.0f} gegen {lpb['beladung_MWh']:,.0f} MWh/a")

c3 = copy.deepcopy(c)
c3["last"]["modus"] = "direkt"
c3["last"]["feldlast_kW"] = list(lp0["P_kW"])
pruefe("modus 'direkt' liefert dasselbe Profil wie 'solar_bedarf'",
       np.allclose(B.lastprofil(c3)["P_kW"], lp0["P_kW"]))

# Last je Sonde muss durch die Sondenzahl des VOLLEN Feldes geteilt werden
cv = copy.deepcopy(c); cv["field"]["symmetrie"] = "voll"
cq = copy.deepcopy(c); cq["field"]["symmetrie"] = "viertel"
pruefe("Last je Sonde ist unabhängig von Voll- oder Viertelmodell",
       np.allclose(B.leistung_je_sonde_W(cv), B.leistung_je_sonde_W(cq)),
       f"{np.abs(B.leistung_je_sonde_W(cv)).max():,.0f} W Spitze")
pruefe("Summe der Sondenlasten ergibt die Feldlast",
       np.allclose(B.leistung_je_sonde_W(cv) * B._n_full(cv) / 1000.0,
                   lp0["P_kW"]))

print("\n" + "=" * 70)
print("2  Geometrie")
print("=" * 70)

pruefe("Viertel enthält ein Viertel der Sonden",
       len(B.sondenpositionen(cq)) * 4 == len(B.sondenpositionen(cv)),
       f"{len(B.sondenpositionen(cq))} von {len(B.sondenpositionen(cv))}")
pruefe("alle Viertelsonden liegen im Quadranten x>0, y>0",
       all(x > 0 and y > 0 for x, y in B.sondenpositionen(cq)))

cu = copy.deepcopy(c); cu["field"].update(n_x_full=5, symmetrie="viertel")
try:
    B.sondenpositionen(cu)
    ok = False
except ValueError as e:
    ok = "gerade" in str(e)
pruefe("ungerade Sondenzahl im Viertelmodell wird abgefangen", ok)

cz = copy.deepcopy(c); cz["field"]["einheitszelle"] = True
x_lo, y_lo, x_hi, y_hi = B.gebiet(cz)
pruefe("Einheitszelle ist genau ein Raster gross",
       abs((x_hi - x_lo) - c["field"]["abstand_m"]) < 1e-9,
       f"{x_hi-x_lo:.1f} m")

d_ein = B.eindringtiefe_m(c)
LX, LY = B.gebiet_halbmasse(c)
fx, fy = B.feldausdehnung(c)
pruefe("Gebietsrand deckt die Eindringtiefe",
       (LX - fx) >= d_ein, f"Rand {LX-fx:.0f} m, Eindringtiefe {d_ein:.0f} m")

print("\n" + "=" * 70)
print("3  Stoffwerte und Nachrechnung")
print("=" * 70)

L = c["layers"][0]; fl = c["fluid"]
e = B.effektive_stoffwerte(L, fl)
n = L["porositaet"]
pruefe("Wärmeleitfähigkeit ist volumengewichtet gemischt",
       abs(e["lambda"] - ((1-n)*L["lambda_s"] + n*fl["lambda"])) < 1e-12)
pruefe("Produkt Dichte mal cp ergibt die volumetrische Kapazität",
       abs(e["dichte"]*e["cp"] - ((1-n)*L["rho_s"]*L["cp_s"]
                                  + n*fl["rho"]*fl["cp"])) < 1e-6)

korr = B.fluid_korrektur(c)
pruefe("Fluidkorrektur ist positiv (Fluid kälter bei Entnahme)",
       korr > 0, f"{korr:.4f} K/(W/m)")

nref = B._n_full(c)
pruefe("Sondenzahl-Formel ist konsistent mit der Linearität",
       abs(B.sondenzahl_fuer_grenze(c, -10.0, 0.0) - 2*nref) < 1e-6,
       f"-10 °C bei {nref} -> {B.sondenzahl_fuer_grenze(c, -10.0, 0.0):.0f} für 0 °C")

print("\n" + "=" * 70)
print("4  Netzautomatik")
print("=" * 70)

ck = klein(); ck["netz"]["ziel_knoten"] = 20_000
mh = B.netzparameter(ck)
pruefe("Knotenbudget wird eingehalten oder gemeldet",
       mh["knoten"] <= ck["netz"]["ziel_knoten"] or mh["vergroebert"] > 1.0,
       f"{mh['knoten']:,} Knoten, Faktor {mh['vergroebert']:.2f}")
pruefe("Wandauflösung bleibt feiner als die Sondenbox",
       mh["groesse_an_sonde_m"] <= ck["borehole"]["box_dx_m"] + 1e-9,
       f"{mh['groesse_an_sonde_m']:.3f} m gegen Box "
       f"{ck['borehole']['box_dx_m']:.2f} m")

print("\n" + "=" * 70)
print("5  Viertel gegen Vollmodell — der eigentliche Test")
print("=" * 70)

tmp = Path(tempfile.mkdtemp(prefix="btes_pruefung_"))
try:
    cq = klein(symmetrie="viertel"); cv = klein(symmetrie="voll")
    print("  rechne Viertel ...")
    dq = lauf(cq, tmp / "viertel")
    print("  rechne Vollfeld ...")
    dv = lauf(cv, tmp / "voll")

    Eq, Vq, Bq = waermeinhalt(cq, dq)
    Ev, Vv, Bv = waermeinhalt(cv, dv)
    pruefe("Gebietsvolumen verhält sich wie 1 zu 4",
           abs(Vv/Vq - 4.0) < 1e-6, f"{Vv/Vq:.6f}")

    L_bh = c["borehole"]["tiefe_fuss_m"] - c["borehole"]["tiefe_kopf_m"]
    soll_q = len(B.sondenpositionen(cq)) * 0.36 * L_bh
    pruefe("vernetztes Sondenvolumen trifft die Vorgabe",
           abs(Bq - soll_q) / soll_q < 1e-6,
           f"{Bq:.3f} gegen {soll_q:.3f} m3")

    abw = abs(Ev/Eq - 4.0) / 4.0 * 100
    pruefe("Wärmeinhalt Voll zu Viertel ist 4,000",
           abw < 0.5, f"{Ev/Eq:.4f}, Abweichung {abw:.3f} %")

    print("\n  Auswertung auf beiden Läufen:")
    for name, cc_, dd in (("Viertel", cq, dq), ("Voll", cv, dv)):
        z, pos, js, kk = B.zeitreihe(cc_, dd)
        pruefe(f"{name}: Sondenzahl aus dem Netz stimmt",
               len(pos) == len(B.sondenpositionen(cc_)),
               f"{len(pos)} Sonden")
        pruefe(f"{name}: Temperaturmatrix hat die richtige Form",
               js.shape == (len(z), len(pos)), f"{js.shape}")
        pruefe(f"{name}: Anfangstemperatur ist der Ausgangszustand",
               abs(js[0].mean() - c["initial"]["T_C"]) < 1e-6,
               f"{js[0].mean():.4f} °C")
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print("\n" + "=" * 70)
n_ok = sum(1 for _, ok, _ in ERGEBNIS if ok)
print(f"ERGEBNIS: {n_ok} von {len(ERGEBNIS)} Prüfungen bestanden")
for name, ok, detail in ERGEBNIS:
    if not ok:
        print(f"  FEHLER: {name}   {detail}")
print("=" * 70)
sys.exit(0 if n_ok == len(ERGEBNIS) else 1)
