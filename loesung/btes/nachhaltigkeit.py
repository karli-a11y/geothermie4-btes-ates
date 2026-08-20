#!/usr/bin/env python3
"""
Nachhaltigkeit — welcher Reduktionsgrad macht das Feld dauerhaft tragfähig?

Der Speicher in der Aufgabenstellung wird netto entleert: der Entnahme steht
nur etwa halb so viel Beladung gegenüber, also kühlt der Untergrund Jahr für
Jahr weiter aus. Mehr Sonden verlangsamen das nur — ein Saisonalspeicher
verschiebt Wärme INNERHALB eines Jahres, ein dauerhaftes Defizit kann er nicht
liefern. Die Stellschraube ist die Bilanz.

Dieses Skript tastet den Reduktionsgrad des Wärmebedarfs ab und beantwortet
für jeden Wert:

  * wie sich Beladung, Entnahme und Jahresbilanz verschieben,
  * wieviel Prozent der heutigen Entnahme noch gedeckt werden,
  * wo die Fluidtemperatur im letzten Betriebsjahr landet,
  * wieviele Sonden nötig wären, um 0 °C bzw. 4 °C zu halten,
  * ob das Feld einschwingt oder weiter wegläuft.

    python nachhaltigkeit.py

-----------------------------------------------------------------------------
WARUM ZWEI LÄUFE FÜR DEN GANZEN BEREICH GENÜGEN
-----------------------------------------------------------------------------
Das Lastprofil ist die Speicherbilanz aus Solarertrag und Bedarf:

    dQ(m) = A_koll * q_sol(m) - g * D_bedarf(m),      g = 1 - Reduktionsgrad

Es hängt damit AFFIN von g ab. Reine Wärmeleitung ist linear, also gilt dieselbe
Zerlegung für die Antwort:

    T(t; g) - T0 = b_S(t) + g * b_D(t)

mit b_S aus einem Lauf mit NUR dem Solareintrag und b_D aus einem Lauf mit NUR
dem Bedarf. Beide werden einmal an der Einheitszelle gerechnet — zusammen gut
zwei Minuten — und danach folgt jeder Reduktionsgrad durch Überlagerung, ohne
eine einzige weitere Simulation.

Die Einheitszelle ist dabei die exakte Lösung für die Mittensonde eines großen
Feldes, also die kritischste Sonde. Für ein endliches Feld ist das Ergebnis
konservativ: bei ausgeglichener Bilanz verliert ein reales Feld über seinen
Umfang Wärme und wird dadurch eher etwas kälter, nicht wärmer (gemessen am
220-Sonden-Feld: 0,16 K).
"""
from __future__ import annotations

import copy
import csv
import glob
import re
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import btes_loesung as B


# ======================================================================
#  EINSTELLUNGEN
# ======================================================================
STUDIE = {
    # Reduktionsgrad des Wärmebedarfs in Prozent.
    #   0  = Bedarf wie in der Aufgabenstellung
    #  23  = etwa dort schließt sich die Jahresbilanz
    #  100 = kein Bedarf mehr
    # Der Solareintrag bleibt in jedem Fall unverändert.
    "reduktionsgrade_prozent": [0, 5, 10, 15, 20, 23, 25, 30, 40],

    # Temperaturgrenzen, gegen die geprüft wird [°C]
    "grenzen_C": [0.0, 4.0],

    # Betriebsjahre. Kürzer heißt schöneres Ergebnis — bei negativer Bilanz
    # wächst die Auskühlung mit der Zeit, die Aussage hängt also am Horizont.
    "jahre": 30,

    # Auflösung der Einheitszelle. 2 Elemente über die Sondenbox reichen für
    # die Aussage; 3 kostet etwa das Doppelte und verschiebt das Ergebnis um
    # gut 0,1 K nach unten.
    "elemente_je_box": 2.0,

    # Zwischenergebnisse der beiden Basisläufe. Sind sie vorhanden, wird
    # nicht neu gerechnet — zum Neurechnen einfach löschen.
    "ordner": "nachhaltigkeit",

    "abbildung": True,
}


# ======================================================================
#  Basisläufe
# ======================================================================
def _zellenlauf(cfg_basis, monatslast_kW, name: str, ordner: Path):
    """Einheitszelle mit einem vorgegebenen Feld-Monatsprofil rechnen.

    Gibt Zeit [a], Boxtemperatur [°C] und den spezifischen Wärmestrom je
    Ausgabezeitpunkt zurück.
    """
    import pyvista as pv
    cache = ordner / f"{name}.csv"
    if cache.exists():
        with open(cache) as fh:
            r = list(csv.DictReader(fh))
        return (np.array([float(x["t_a"]) for x in r]),
                np.array([float(x["T_box_C"]) for x in r]),
                np.array([float(x["q_W_m"]) for x in r]))

    c = copy.deepcopy(cfg_basis)
    c["ablauf"] = {k: False for k in c["ablauf"]}
    c["field"]["einheitszelle"] = True
    c["netz"]["elemente_je_box"] = STUDIE["elemente_je_box"]
    c["zeit"]["jahre"] = STUDIE["jahre"]
    c["last"]["modus"] = "direkt"
    c["last"]["feldlast_kW"] = list(monatslast_kW)
    d = ordner / name
    d.mkdir(parents=True, exist_ok=True)
    c["ausgabe"] = dict(c["ausgabe"]); c["ausgabe"]["ordner"] = str(d)

    t0 = time.time()
    B.netz_bauen(c, d)
    nf = B.netz_wandeln(c, d / f"{c['ausgabe']['prefix']}.msh", d)
    rc = B.ogs_starten(B.prj_bauen(c, d, nf, B.lastkurve(c)))
    if rc != 0:
        raise SystemExit(f"OGS endete mit {rc} in {d}")
    print(f"    {name}: {time.time()-t0:.0f} s")

    dom = pv.read(str(d / f"{c['ausgabe']['prefix']}_domain.vtu"))
    mid = np.asarray(dom.cell_data["MaterialIDs"])
    vol = np.abs(np.asarray(dom.compute_cell_sizes(
        length=False, area=False, volume=True).cell_data["Volume"]))
    sel = mid >= len(c["layers"])
    w = vol[sel] / vol[sel].sum()
    P = B.leistung_je_sonde_W(c)
    L = c["borehole"]["tiefe_fuss_m"] - c["borehole"]["tiefe_kopf_m"]

    t, T, q = [], [], []
    for f in sorted(glob.glob(str(d / f"{c['ausgabe']['prefix']}_ts_*.vtu")),
                    key=lambda x: float(re.search(r"_t_([\d.]+)\.vtu$",
                                                  x).group(1))):
        td = float(re.search(r"_t_([\d.]+)\.vtu$", f).group(1)) / B.DAY
        m = pv.read(f)
        Tc = np.asarray(m.point_data_to_cell_data().cell_data["T"])
        t.append(td / B.YEAR_D)
        T.append(float((Tc[sel] * w).sum()) - 273.15)
        q.append(P[(int(round(td / B.MONAT_D)) - 1) % 12] / L)
    t, T, q = np.array(t), np.array(T), np.array(q)
    with open(cache, "w", newline="") as fh:
        wr = csv.writer(fh); wr.writerow(["t_a", "T_box_C", "q_W_m"])
        wr.writerows(zip(t, T, q))
    return t, T, q


def basislaeufe(cfg):
    """Solareintrag allein und Bedarf allein — die Bausteine für alles Weitere."""
    L = cfg["last"]
    if L["modus"] != "solar_bedarf":
        raise SystemExit("Diese Studie braucht last.modus = 'solar_bedarf'.")
    q = np.asarray(L["solarertrag_kWh_m2_monat"], float)
    prof = np.asarray(L["bedarfsprofil"], float); prof = prof / prof.sum()
    E_sol = L["kollektorflaeche_m2"] * q                       # kWh/Monat
    E_bed = L["waermebedarf_MWh_a"] * 1e3 * prof               # kWh/Monat, g=1
    P_sol = E_sol / (B.MONAT_D * 24.0)                         # kW
    P_bed = E_bed / (B.MONAT_D * 24.0)

    ordner = HERE / STUDIE["ordner"]
    ordner.mkdir(parents=True, exist_ok=True)
    print("Basisläufe an der Einheitszelle "
          f"({STUDIE['jahre']} Jahre, {STUDIE['elemente_je_box']:.0f} "
          f"Elemente je Box):")
    t_s, b_s, q_s = _zellenlauf(cfg, +P_sol, "nur_solar", ordner)
    t_d, b_d, q_d = _zellenlauf(cfg, -P_bed, "nur_bedarf", ordner)
    if not np.allclose(t_s, t_d):
        raise SystemExit("Basisläufe haben unterschiedliche Ausgabezeiten.")
    T0 = cfg["initial"]["T_C"]
    return dict(t=t_s, bS=b_s - T0, bD=b_d - T0, qS=q_s, qD=q_d,
                E_sol=E_sol, E_bed=E_bed, T0=T0)


# ======================================================================
#  Überlagerung
# ======================================================================
def kennzahlen(cfg, bas, reduktion_prozent):
    """Alle Kennzahlen für einen Reduktionsgrad — ohne neue Simulation."""
    g = 1.0 - reduktion_prozent / 100.0
    dQ = bas["E_sol"] - g * bas["E_bed"]              # kWh/Monat
    E = dQ / 1e3                                       # MWh/Monat
    lade, entl = float(E[E > 0].sum()), float(-E[E < 0].sum())

    korr = B.fluid_korrektur(cfg)
    fluid = bas["T0"] + bas["bS"] + g * bas["bD"] + (bas["qS"] + g * bas["qD"]) * korr
    t = bas["t"]
    letzte = t > t.max() - 1.0
    vorletzte_dekade = (t > t.max() - 11.0) & (t <= t.max() - 10.0)
    tmin = float(fluid[letzte].min())
    drift = (tmin - float(fluid[vorletzte_dekade].min())
             if vorletzte_dekade.any() else float("nan"))

    L = cfg["borehole"]["tiefe_fuss_m"] - cfg["borehole"]["tiefe_kopf_m"]
    P = dQ / (B.MONAT_D * 24.0) * 1000.0 / B._n_full(cfg)      # W je Sonde
    return dict(reduktion=reduktion_prozent, g=g, beladung=lade, entnahme=entl,
                bilanz=lade - entl,
                anteil=100.0 * entl / lade if lade > 0 else float("inf"),
                tmin=tmin, drift=drift, spez=float(np.abs(P).max()) / L,
                fluid=fluid)


def _grenzsuche(cfg, bas, ziel_bilanz=0.0):
    """Reduktionsgrad, bei dem sich die Jahresbilanz genau schließt."""
    g = float(bas["E_sol"].sum()) / float(bas["E_bed"].sum())
    return (1.0 - g) * 100.0


# ======================================================================
def main() -> int:
    cfg = B.CONFIG
    bas = basislaeufe(cfg)
    n_ist = B._n_full(cfg)
    r_bal = _grenzsuche(cfg, bas)
    ref = kennzahlen(cfg, bas, 0.0)

    print(f"\nAusgangslage: {n_ist} Sonden, Entnahme "
          f"{ref['entnahme']:,.0f} MWh/a, Bilanz {ref['bilanz']:+,.0f} MWh/a")
    print(f"Bilanz schließt sich bei Reduktionsgrad {r_bal:.1f} %\n")

    grade = sorted(set(STUDIE["reduktionsgrade_prozent"] + [round(r_bal, 1)]))
    zeilen = [kennzahlen(cfg, bas, r) for r in grade]

    kopf = (f"{'Red.':>5} {'Bedarf':>8} {'Belad.':>8} {'Entn.':>8} "
            f"{'Entn/Inj':>9} {'Bilanz':>8} {'W/m':>6} {'Fluid':>8} "
            f"{'K/10a':>7}")
    for gz in STUDIE["grenzen_C"]:
        kopf += f" {'N(' + f'{gz:.0f}' + '°C)':>9}"
    print(kopf)
    print(f"{'%':>5} {'MWh/a':>8} {'MWh/a':>8} {'MWh/a':>8} {'%':>9} "
          f"{'MWh/a':>8} {'':>6} {'°C':>8} {'':>7}"
          + "".join(f" {'':>9}" for _ in STUDIE["grenzen_C"]))
    print("-" * len(kopf))
    for k in zeilen:
        bed = cfg["last"]["waermebedarf_MWh_a"] * k["g"]
        zeile = (f"{k['reduktion']:>5.1f} {bed:>8,.0f} {k['beladung']:>8,.0f} "
                 f"{k['entnahme']:>8,.0f} {k['anteil']:>9.1f} "
                 f"{k['bilanz']:>+8,.0f} {k['spez']:>6.1f} {k['tmin']:>8.2f} "
                 f"{k['drift']:>+7.2f}")
        for gz in STUDIE["grenzen_C"]:
            n = B.sondenzahl_fuer_grenze(cfg, k["tmin"], gz)
            zeile += f" {n:>9,.0f}" if 0 < n < 1e6 else f" {'—':>9}"
        marke = "  <-- Bilanz null" if abs(k["reduktion"] - r_bal) < 0.06 else ""
        print(zeile + marke)

    print(f"\nN(...) = nötige Sondenzahl für diese Grenze, vorhanden {n_ist}.")
    print("Drift = Änderung des Jahresminimums über das letzte Jahrzehnt.")
    print("        Nahe null heißt eingeschwungen; deutlich negativ heißt,")
    print("        das Feld läuft weiter weg und der Wert im Jahr "
          f"{STUDIE['jahre']} ist nur eine Momentaufnahme.")

    with open(HERE / "nachhaltigkeit.csv", "w", newline="") as fh:
        felder = ["reduktion", "g", "beladung", "entnahme", "bilanz", "anteil",
                  "tmin", "drift", "spez"]
        wr = csv.DictWriter(fh, fieldnames=felder, extrasaction="ignore")
        wr.writeheader(); wr.writerows(zeilen)
    print(f"\n-> {HERE / 'nachhaltigkeit.csv'}")

    if STUDIE["abbildung"]:
        _abbildung(cfg, bas, zeilen, r_bal, n_ist)
    return 0


def _abbildung(cfg, bas, zeilen, r_bal, n_ist):
    plt = B._plt(); F = B._fmt()
    fein = [kennzahlen(cfg, bas, r) for r in np.linspace(0, 45, 120)]
    r = np.array([k["reduktion"] for k in fein])
    tm = np.array([k["tmin"] for k in fein])
    bil = np.array([k["bilanz"] for k in fein])
    ant = np.array([k["anteil"] for k in fein])

    fig, ax = plt.subplots(1, 3, figsize=(14.0, 4.5))

    a = ax[0]
    a.axhline(0, color=B.AXIS, lw=1.0)
    a.plot(r, bil, color=B.C1)
    a.axvline(r_bal, color=B.C3, lw=1.2, ls="--")
    a.annotate(f"Bilanz null bei\n{B.de(r_bal,1)} %", (r_bal, bil.min() * 0.55),
               color=B.C3, fontsize=8.5, fontweight="bold", ha="left",
               xytext=(6, 0), textcoords="offset points")
    a.set_xlabel("Reduktionsgrad des Wärmebedarfs [%]")
    a.set_ylabel("Jahresbilanz [MWh/a]")
    a.set_title("Energiebilanz des Speichers")
    a.xaxis.set_major_formatter(F); a.yaxis.set_major_formatter(F)
    a.grid(alpha=0.9); a.set_axisbelow(True)

    b = ax[1]
    for gz, col in zip(STUDIE["grenzen_C"], (B.C4, B.C2)):
        b.axhline(gz, color=col, lw=1.2, ls="--")
        b.annotate(f"{B.de(gz,0)} °C", (0.5, gz), color=col, fontsize=8.5,
                   fontweight="bold", va="bottom",
                   xytext=(0, 3), textcoords="offset points")
    b.plot(r, tm, color=B.C1)
    b.axvline(r_bal, color=B.C3, lw=1.2, ls="--")
    b.set_xlabel("Reduktionsgrad des Wärmebedarfs [%]")
    b.set_ylabel(f"Fluidminimum im Jahr {STUDIE['jahre']} [°C]")
    b.set_title(f"Mit den vorhandenen {n_ist} Sonden")
    b.set_ylim(max(-40, tm.min() * 1.05), max(8, tm.max() * 1.1))
    b.xaxis.set_major_formatter(F); b.yaxis.set_major_formatter(F)
    b.grid(alpha=0.9); b.set_axisbelow(True)

    c = ax[2]
    c.axhline(100, color=B.AXIS, lw=1.0, ls="--")
    c.annotate("Entnahme = Injektion", (0.5, 100), color=B.INK2, fontsize=8.5,
               va="bottom", xytext=(0, 3), textcoords="offset points")
    c.plot(r, ant, color=B.C2)
    c.axvline(r_bal, color=B.C3, lw=1.2, ls="--")
    c.set_xlabel("Reduktionsgrad des Wärmebedarfs [%]")
    c.set_ylabel("Entnahme in % der Beladung")
    c.set_title("Wie stark der Speicher überzogen wird")
    c.xaxis.set_major_formatter(F); c.yaxis.set_major_formatter(F)
    c.grid(alpha=0.9); c.set_axisbelow(True)

    fig.suptitle(f"Nachhaltigkeit — der Speicher trägt dauerhaft ab einem "
                 f"Reduktionsgrad von {B.de(r_bal,1)} %",
                 x=0.006, ha="left", fontsize=11.5, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    ziel = B._figdir(cfg) / "7_nachhaltigkeit.png"
    fig.savefig(ziel, dpi=150); plt.close(fig)
    print(f"-> {ziel}")


if __name__ == "__main__":
    sys.exit(main())
