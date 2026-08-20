#!/usr/bin/env python3
"""
BTES-Sondenfeld — Modell, Auswertung und Abbildungen in einer Datei.

Rechnet ein Erdsonden-Wärmespeicherfeld als reine Wärmeleitung, wertet den Lauf
aus und erzeugt die Abbildungen. Gesteuert wird alles über CONFIG ganz oben;
es gibt keine Kommandozeilenschalter.

    python btes_loesung.py

Was gerechnet wird, steht im Block "ablauf": Bilanz, Netz, OGS, Auswertung,
Abbildungen lassen sich einzeln an- und abschalten. Für einen ersten Blick
genügt "bilanz" allein — das dauert Sekunden und zeigt Energiebilanz,
spezifische Entzugsrate und die Modellgröße, bevor irgendetwas stundenlang
rechnet.

-----------------------------------------------------------------------------
DREI VEREINFACHUNGEN
-----------------------------------------------------------------------------
1) Prozess T statt HT. Das Druckfeld ist bei diesem Aufbau identisch null —
   kein Auftrieb, keine Schwerkraft, p = 0 oben und unten, seitlich kein Rand.
   Die Druckgleichung ist homogen, ihre einzige Lösung p = 0. Damit verschwindet
   die Darcy-Geschwindigkeit, der Advektionsterm entfällt, übrig bleibt reine
   Wärmeleitung. Halbiert die Unbekannten.
   Auch mit einem regionalen Gradienten i = 2e-3 bliebe die Wärmeleitung
   dominant: in der durchlässigsten Schicht (k ~ 5e-14 m2) ergibt sich
   v ~ 1e-9 m/s, also 3 cm im Jahr, und eine Péclet-Zahl von 0,03. Relevant
   würde Strömung erst bei etwa dem Achtzigfachen der Durchlässigkeit.

2) Viertelmodell (optional). Feld, Schichten, Last und Randbedingungen sind
   spiegelsymmetrisch zu x = 0 und y = 0. Auf einer Spiegelebene verschwindet
   der Wärmestrom senkrecht dazu — und eine Fläche OHNE Randbedingung ist in
   OGS genau das. Die Schnittflächen bekommen deshalb bewusst keine
   Randbedingung; die Symmetriebedingung ist damit exakt erfüllt. Faktor 4.
   Voraussetzung: gerade Sondenzahl je Richtung, sonst läge eine Sonde auf der
   Schnittebene. Nachgewiesen: Wärmeinhalt Voll zu Viertel = 4,0000 bei
   0,001 % Abweichung.

3) Prismennetz. Das Temperaturfeld um eine Sonde ist quasi eindimensional
   radial; in der Tiefe passiert außer an Kopf und Fuß fast nichts. Ein
   2D-Grundriss wird in z extrudiert: fein in der Ebene, grob in der Tiefe,
   mit konformen Schichtgrenzen.

Zusammen etwa Faktor 130 je Zeitschritt gegenüber dem HT-Vollmodell.

-----------------------------------------------------------------------------
LAST — die häufigste Fehlerquelle
-----------------------------------------------------------------------------
Angegeben wird die GESAMTLAST DES FELDES, nie die Last je Sonde. Das Skript
teilt intern durch die Sondenzahl des VOLLEN Feldes und prägt das Ergebnis
jeder modellierten Sonde auf. Wer selbst durch die Sondenzahl des Viertels
teilt, rechnet mit der vierfachen Last.

-----------------------------------------------------------------------------
REFERENZFÄLLE ZUM NACHRECHNEN
-----------------------------------------------------------------------------
Vier geprüfte Kombinationen aus Sondenzahl und Last. Zum Nachrechnen genügt es,
drei Zeilen in CONFIG zu setzen und field.einheitszelle = True zu lassen — dann
ist jeder Fall in gut einer Minute durch. Angegeben ist das tiefste
Fluidminimum im dreißigsten Betriebsjahr.

  Fall  n_x  n_y  Reduktion   W/m   Fluid J30   Ergebnis
  ----  ---  ---  ---------  -----  ---------  ------------------------------
   A     22   10     0 %      36.0    -55.18   friert ein
   B     22   10    23 %      25.9     +0.42   hält 0 °C, keine Reserve
   C     26   14    23 %      15.7     +4.21   hält 4 °C
   D     22   10    25 %      25.0     +5.27   hält 4 °C, gleiche Sondenzahl

Fall A ist die Aufgabenstellung wie geliefert: der Entnahme steht nur halb so
viel Beladung gegenüber, das Feld kühlt Jahr für Jahr weiter aus. Keine
Sondenzahl repariert das — ein Saisonalspeicher verschiebt Wärme INNERHALB
eines Jahres, ein dauerhaftes Defizit kann er nicht liefern.

Fall B schließt die Jahresbilanz allein über den Bedarf. Die 23 % ergeben sich
aus Solarertrag geteilt durch Bedarf; das Skript rechnet den Wert selbst aus und
meldet ihn bei jedem Lauf. Damit genügen die vorhandenen 220 Sonden — aber nur
knapp, mit 0,42 K Abstand zur Vereisung.

Fall C und D erreichen beide die Fördergrenze von 4 °C, auf zwei verschiedenen
Wegen: C über mehr Bohrungen bei gleichem Bedarf, D über zwei Prozentpunkte
mehr Reduktion bei gleicher Sondenzahl. **Zwei Prozentpunkte Reduktion ersetzen
144 Bohrungen** — das ist die wirtschaftlich entscheidende Zahl.

Der Vergleich zeigt außerdem den Zielkonflikt bei der Auslastung. B und D
liegen mit 25,9 bzw. 25,0 W/m mitten im Literaturbereich für Erdwärmesonden
(20 bis 70 W/m). C nutzt die Sonden mit 15,7 W/m nur noch schwach aus — man
bezahlt Bohrungen, die kaum belastet werden.

-----------------------------------------------------------------------------
RECHENZEITEN (üblicher Arbeitsplatzrechner)
-----------------------------------------------------------------------------
    Einheitszelle, 30 Jahre                    Sekunden bis Minuten
    Viertelfeld 220 Sonden, 30 Jahre           etwa 6 h
    Vollfeld 220 Sonden, 30 Jahre              etwa 80 h

Für Parameterstudien immer zuerst die Einheitszelle: sie ist die exakte Lösung
für die Mittensonde eines großen Feldes. Und weil reine Wärmeleitung linear in
der Last ist, folgt aus EINEM Lauf die nötige Sondenzahl für jede
Temperaturgrenze in geschlossener Form.
"""
from __future__ import annotations

import csv
import glob
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

HERE = Path(__file__).resolve().parent
DAY, YEAR_D = 86400.0, 365.25
MONAT_D = YEAR_D / 12.0
MONATE = ["Jan", "Feb", "Mrz", "Apr", "Mai", "Jun",
          "Jul", "Aug", "Sep", "Okt", "Nov", "Dez"]


# ======================================================================
#  CONFIG
# ======================================================================
CONFIG: dict = {

    # ==================================================================
    #  WAS GERECHNET WIRD
    # ==================================================================
    "ablauf": {
        "bilanz":        True,   # Lastprofil, Energiebilanz, Modellgröße
        "netz":          True,   # Netz und Projektdatei erzeugen
        "ogs":           True,   # Simulation starten
        "auswertung":    True,   # Vereisungsprüfung aus der Ausgabe
        "abbildungen":   True,   # Geometrie, Untergrund, Lastprofil
        "beispiellauf":  False,  # zwei Einheitszellenläufe für Abb. 4 und 5
        "pdf":           False,  # LOESUNG.md nach PDF (braucht pandoc+Chrome)
    },

    # ==================================================================
    #  DAS WICHTIGSTE ZUERST — Last, Feld, Zeitraum
    #
    #  Darunter bleibt alles einstellbar. Untergrund und Prozess sind
    #  vorgegeben: gerechnet wird reine Wärmeleitung. Gebiet und Netz leiten
    #  sich selbst aus dem Feld und der Laufzeit ab.
    # ==================================================================

    # ------------------------------------------------------------------
    # LAST
    # ------------------------------------------------------------------
    # modus "solar_bedarf": das Feld sieht die Speicherbilanz
    #     dQ(m) = kollektorflaeche * solarertrag(m) - g * bedarf(m)
    #     mit g = 1 - reduktionsgrad_prozent / 100
    # modus "direkt": feldlast_kW wird unverändert benutzt (12 Werte, Feld
    #     GESAMT in kW, positiv = Beladung).
    "last": {
        "modus": "solar_bedarf",
        "kollektorflaeche_m2": 6317.0,
        "solarertrag_kWh_m2_monat": [31.5, 48.1, 74.8, 103.0, 117.7, 116.2,
                                     122.2, 115.4, 88.7, 65.6, 33.0, 23.3],
        "waermebedarf_MWh_a": 7707.0,
        "bedarfsprofil": [0.1437, 0.1368, 0.1150, 0.0813, 0.0616, 0.0440,
                          0.0310, 0.0288, 0.0570, 0.0949, 0.0968, 0.1092],

        # ---- DIE ZENTRALE STELLSCHRAUBE -------------------------------
        # Um wieviel Prozent der Wärmebedarf gesenkt wird. Der Solareintrag
        # bleibt dabei unverändert. Ein Wert je Lauf.
        #    0 %  Bedarf wie in der Aufgabenstellung — das Feld kühlt aus
        #   23 %  Jahresbilanz schließt sich, 220 Sonden halten 0 °C
        #   25 %  220 Sonden halten auch die Fördergrenze von 4 °C
        # Das Skript rechnet den Wert für eine ausgeglichene Bilanz selbst aus
        # und meldet ihn bei jedem Lauf.
        "reduktionsgrad_prozent": 0.0,
        "feldlast_kW": None,
        "rampe_tage": 7.0,
    },

    # ------------------------------------------------------------------
    # SONDENFELD
    # ------------------------------------------------------------------
    # Referenzfälle siehe Kopf: A = 22/10 mit 0 % (friert ein),
    # B = 22/10 mit 23 % (hält 0 °C), C = 26/14 mit 23 %,
    # D = 22/10 mit 25 % (hält 4 °C bei gleicher Sondenzahl).
    "field": {
        "n_x_full": 22,
        "n_y_full": 10,
        "abstand_m": 6.0,
        # "voll" oder "viertel" — beide liefern dieselbe Lösung, das Viertel
        # kostet ein Viertel der Unbekannten und rund ein Sechstel der Zeit.
        "symmetrie": "voll",
        # True ersetzt das Feld durch EINE Zelle mit adiabaten Rändern: die
        # exakte Lösung für die Mittensonde eines großen Feldes. Für
        # Parameterstudien. Die Feldlast wird weiterhin durch n_x*n_y geteilt.
        "einheitszelle": False,
    },

    # ------------------------------------------------------------------
    # ZEITRAUM
    # ------------------------------------------------------------------
    "zeit": {
        "jahre": 30,
        "schritte_je_monat": 4,
        "ausgabe_je_n_schritte": 4,
    },

    # ==================================================================
    #  AB HIER: seltener zu ändern — aber alles bleibt offen
    # ==================================================================

    "borehole": {
        "tiefe_kopf_m":  2.0,
        "tiefe_fuss_m": 159.0,
        "box_dx_m": 0.6,
        "box_dy_m": 0.6,
        # nur für die Nachrechnung der Fluidtemperatur, nicht für die Lösung
        "r_bohrloch_m": 0.075,
        "R_b_Km_per_W": 0.10,
    },

    # Gebiet: mindestens die Eindringtiefe sqrt(4*alpha*t) jenseits der
    # äußersten Sonde, sonst sperrt der adiabate Rand die Störung ein.
    "domain": {
        "automatisch": True,
        "rand_faktor": 2.0,
        "rand_min_m": 40.0,
        "size_x_m": None,
        "size_y_m": None,
        "z_basis_m": 0.0,
    },

    # Netz: empfindlich ist ausschließlich die Ebene an der Sondenwand.
    # dz von 4 m auf 1 m verschiebt das Ergebnis um 14 mK, die Ebene von
    # 0,40 m auf 0,20 m um 73 mK. Angestrebt wird box_dx/3; reißt das
    # Ergebnis das Knotenbudget, vergröbert das Skript und sagt es.
    # Obergrenze für die Modellgröße. Gemessen am Vollfeld mit 220 Sonden
    # (Laufzeit skaliert mit n^1,86, Fehler aus der Netzstudie oben):
    #
    #   Budget   Knoten   Wand [m]    Fehler    1 Jahr   30 Jahre
    #    1,5 M    1,50 M    0,478    > 416 mK     8 min      4,1 h
    #    2,5 M    2,51 M    0,358      378 mK    22 min     10,5 h
    #    4,0 M    4,00 M    0,279      292 mK    51 min     25,1 h
    #    9,0 M    7,35 M    0,200      163 mK     2,6 h     77,8 h
    #
    # Der Fehler gilt für den Fall wie gegeben mit 66 K Auslenkung; bei
    # ausgeglichener Bilanz (rund 10 K) ist er etwa ein Sechstel davon.
    # ALLE Netze rechnen zu warm — ein grobes Netz verschmiert die
    # Volumenquelle und macht die Sonde weniger extrem. Für einen knappen
    # Fall taugt ein grobes Netz deshalb nie als Nachweis.
    #
    # Unterhalb von 0,40 m Wandauflösung liegt keine Messung vor; die
    # 1,5-Mio.-Zeile ist nach unten offen und nur für einen ersten Blick.
    "netz": {
        "automatisch": True,
        "ziel_knoten": 2_500_000,
        "elemente_je_box": 3.0,
        "dz_je_sonde": 40,
        "groesse_an_sonde_m": None, "ring_radius_m": None,
        "groesse_im_feld_m": None, "feld_radius_m": None,
        "groesse_fern_m": None, "fern_radius_m": None,
        "dz_sondenbereich_m": None, "dz_ausserhalb_m": None,
    },

    # Untergrund, von OBEN nach UNTEN. Für einen homogenen Boden genügt ein
    # einziger Eintrag; die Sonde darf beliebig viele Schichten durchstoßen.
    "layers": [
        {"name": "Unterer Buntsandstein (suC)", "maechtigkeit_m": 6.0,
         "porositaet": 0.1542, "rho_s": 2630.0, "cp_s": 608.4, "lambda_s": 2.71},
        {"name": "Tonstein Obere Broeckelschiefer (Z8B2T)", "maechtigkeit_m": 9.0,
         "porositaet": 0.0485, "rho_s": 2710.0, "cp_s": 741.7, "lambda_s": 2.48},
        {"name": "Sandstein Obere Broeckelschiefer (Z8B2S)", "maechtigkeit_m": 9.0,
         "porositaet": 0.1484, "rho_s": 2700.0, "cp_s": 585.2, "lambda_s": 1.68},
        {"name": "Tonstein Untere Broeckelschiefer (Z7B1T)", "maechtigkeit_m": 18.0,
         "porositaet": 0.0297, "rho_s": 2700.0, "cp_s": 755.6, "lambda_s": 2.64},
        {"name": "Sandstein Untere Broeckelschiefer (Z7B1S)", "maechtigkeit_m": 10.0,
         "porositaet": 0.1925, "rho_s": 2720.0, "cp_s": 584.6, "lambda_s": 1.57},
        {"name": "Friesland-Ton (T6)", "maechtigkeit_m": 5.0,
         "porositaet": 0.0928, "rho_s": 2770.0, "cp_s": 707.6, "lambda_s": 2.40},
        {"name": "Friesland-Sandstein (S6)", "maechtigkeit_m": 5.0,
         "porositaet": 0.2692, "rho_s": 2680.0, "cp_s": 496.3, "lambda_s": 1.01},
        {"name": "Unterer Ohre-Ton (z5T)", "maechtigkeit_m": 3.0,
         "porositaet": 0.1227, "rho_s": 2730.0, "cp_s": 648.4, "lambda_s": 1.88},
        {"name": "Ohre-Sandstein (S5)", "maechtigkeit_m": 3.0,
         "porositaet": 0.0983, "rho_s": 2680.0, "cp_s": 619.4, "lambda_s": 2.55},
        {"name": "Unterer Aller-Ton (z4T)", "maechtigkeit_m": 3.0,
         "porositaet": 0.2046, "rho_s": 2900.0, "cp_s": 689.7, "lambda_s": 1.83},
        {"name": "Aller-Sandstein (S4)", "maechtigkeit_m": 7.0,
         "porositaet": 0.2009, "rho_s": 2650.0, "cp_s": 558.5, "lambda_s": 2.54},
        {"name": "Oberer Leine-Ton (T3r)", "maechtigkeit_m": 4.0,
         "porositaet": 0.1500, "rho_s": 2500.0, "cp_s": 880.0, "lambda_s": 2.20},
        {"name": "Plattendolomit (Ca3)", "maechtigkeit_m": 38.0,
         "porositaet": 0.1191, "rho_s": 2740.0, "cp_s": 850.4, "lambda_s": 3.32},
        {"name": "Unterer Leine-Ton (T3)", "maechtigkeit_m": 3.0,
         "porositaet": 0.0772, "rho_s": 2850.0, "cp_s": 663.2, "lambda_s": 1.87},
        {"name": "Oberer Strassfurt-Ton (T2r)", "maechtigkeit_m": 20.0,
         "porositaet": 0.1500, "rho_s": 2500.0, "cp_s": 880.0, "lambda_s": 2.20},
        {"name": "Strassfurt-Sulfat (A2T)", "maechtigkeit_m": 10.0,
         "porositaet": 0.1455, "rho_s": 2680.0, "cp_s": 750.0, "lambda_s": 1.05},
        {"name": "Hauptdolomit (Ca2)", "maechtigkeit_m": 4.0,
         "porositaet": 0.2158, "rho_s": 2760.0, "cp_s": 702.9, "lambda_s": 2.40},
        {"name": "Oberer Werra-Ton (T1r)", "maechtigkeit_m": 21.0,
         "porositaet": 0.1500, "rho_s": 2500.0, "cp_s": 880.0, "lambda_s": 2.20},
        {"name": "Oberer Werra-Anhydrit (A1r)", "maechtigkeit_m": 21.0,
         "porositaet": 0.2553, "rho_s": 2690.0, "cp_s": 446.1, "lambda_s": 0.56},
        {"name": "Unterer Werra-Anhydrit (A1)", "maechtigkeit_m": 70.0,
         "porositaet": 0.0179, "rho_s": 2950.0, "cp_s": 633.9, "lambda_s": 4.48},
        {"name": "Anhydritknotenschiefer (A1Ca)", "maechtigkeit_m": 3.0,
         "porositaet": 0.2657, "rho_s": 2670.0, "cp_s": 644.2, "lambda_s": 0.75},
        {"name": "Zechsteinkalk (Ca1)", "maechtigkeit_m": 3.2,
         "porositaet": 0.0114, "rho_s": 2690.0, "cp_s": 806.7, "lambda_s": 2.56},
    ],

    # Porenwasser. Geht nur über die Mischung in die Stoffwerte ein.
    "fluid": {"rho": 1000.0, "cp": 4180.0, "lambda": 0.6},

    "initial": {"T_C": 10.0, "geothermischer_gradient_K_m": 0.0},

    "ausgabe": {"prefix": "btes", "ordner": "out", "figures": "figures"},

    "loeser": {
        # Reine Wärmeleitung -> Systemmatrix symmetrisch positiv definit.
        "typ": "CG", "vorkonditionierer": "DIAGONAL",
        "toleranz": 1.0e-10, "max_iter": 20000,
        "nichtlinear_max_iter": 20, "reltol_T": 1.0e-2,
    },
}

GRENZEN = [(0.0, "Vereisung 0 °C"), (4.0, "Fördergrenze 4 °C")]


# ======================================================================
#  Last und Stoffwerte
# ======================================================================
def lastprofil(cfg: dict) -> dict:
    """Monatliches Lastprofil des FELDES [kW]; positiv = Beladung."""
    L = cfg["last"]
    if L["modus"] == "direkt":
        P = np.asarray(L["feldlast_kW"], float)
        if P.size != 12:
            raise ValueError("last.feldlast_kW braucht genau 12 Werte.")
        E_sol = E_bed = None
    else:
        q = np.asarray(L["solarertrag_kWh_m2_monat"], float)
        prof = np.asarray(L["bedarfsprofil"], float)
        if q.size != 12 or prof.size != 12:
            raise ValueError("Solarertrag und Bedarfsprofil brauchen 12 Werte.")
        prof = prof / prof.sum()
        g = 1.0 - float(L.get("reduktionsgrad_prozent", 0.0)) / 100.0
        E_sol = L["kollektorflaeche_m2"] * q
        E_bed = g * L["waermebedarf_MWh_a"] * 1e3 * prof
        P = (E_sol - E_bed) / (MONAT_D * 24.0)

    E = P * MONAT_D * 24.0 / 1e3
    lade, entl = float(E[E > 0].sum()), float(-E[E < 0].sum())
    out = {"P_kW": P, "beladung_MWh": lade, "entnahme_MWh": entl,
           "bilanz_MWh": lade - entl}
    if E_sol is not None:
        out["solar_MWh"] = float(E_sol.sum()) / 1e3
        out["bedarf_MWh"] = float(E_bed.sum()) / 1e3
        # Reduktionsgrad, bei dem sich die Jahresbilanz genau schließt
        g_null = float(E_sol.sum()) / (L["waermebedarf_MWh_a"] * 1e3)
        out["faktor_bilanz_null"] = g_null
        out["reduktion_bilanz_null"] = (1.0 - g_null) * 100.0
    return out


def _n_full(cfg) -> int:
    return int(cfg["field"]["n_x_full"]) * int(cfg["field"]["n_y_full"])


def leistung_je_sonde_W(cfg) -> np.ndarray:
    """Feldlast durch die Sondenzahl des VOLLEN Feldes — nie durch die des
    Viertels. Das ist die häufigste Fehlerquelle beim Symmetriemodell."""
    return lastprofil(cfg)["P_kW"] * 1000.0 / _n_full(cfg)


def effektive_stoffwerte(mat: dict, fluid: dict) -> dict:
    """Korn und Porenwasser mischen — genau wie der HT-Prozess es intern tut.

    lambda_eff = (1-n)*lambda_s + n*lambda_f
    (rho c)_eff = (1-n)*rho_s*cp_s + n*rho_f*cp_f
    In die Gleichung geht nur das Produkt Dichte mal cp ein; die Aufteilung
    ist frei und hier so gewählt, dass beide Zahlen lesbar bleiben.
    """
    n = float(mat["porositaet"])
    rho_eff = (1 - n) * mat["rho_s"] + n * fluid["rho"]
    rhoc = (1 - n) * mat["rho_s"] * mat["cp_s"] + n * fluid["rho"] * fluid["cp"]
    return {"dichte": rho_eff, "cp": rhoc / rho_eff,
            "lambda": (1 - n) * mat["lambda_s"] + n * fluid["lambda"],
            "rho_c": rhoc}


def fluid_korrektur(cfg) -> float:
    """Faktor [K/(W/m)] von der Boxtemperatur zur Fluidtemperatur.

        T_Fluid = T_Box + q' * [ ln(r_box/r_b)/(2*pi*lambda) + R_b ]

    Die 0,6-m-Box ist kein Bohrloch: sie unterschätzt den Abfall zur Wand, und
    der Bohrlochwiderstand fehlt ganz. Bei Entnahme (q' < 0) ist das Fluid
    also kälter als die Box — hier um rund 5 K. Nachrechnung, keine
    Modellgröße; r_b und R_b sind Annahmen und gehen linear ein.
    """
    bh = cfg["borehole"]
    r_box = np.sqrt(bh["box_dx_m"] * bh["box_dy_m"] / np.pi)
    mitte = 0.5 * (bh["tiefe_kopf_m"] + bh["tiefe_fuss_m"])
    d, lam = 0.0, cfg["layers"][0]["lambda_s"]
    for L in cfg["layers"]:
        d += L["maechtigkeit_m"]
        if d >= mitte:
            lam = L["lambda_s"]
            break
    return (np.log(r_box / bh["r_bohrloch_m"]) / (2 * np.pi * lam)
            + bh["R_b_Km_per_W"])


def sondenzahl_fuer_grenze(cfg, T_min_C, grenze_C) -> float:
    """Nötige Sondenzahl, damit die Fluidtemperatur grenze_C gerade hält.

    Reine Wärmeleitung ist linear in der Last; bei fester Feldlast ist die
    Last je Sonde proportional zu 1/N, also
        T(N) - T0 = (T(N_ref) - T0) * N_ref / N
    """
    T0 = cfg["initial"]["T_C"]
    if grenze_C >= T0 or T_min_C >= grenze_C:
        return float(_n_full(cfg))
    return _n_full(cfg) * (T_min_C - T0) / (grenze_C - T0)


# ======================================================================
#  Geometrie
# ======================================================================
def _viertel(cfg) -> bool:
    return cfg["field"].get("symmetrie", "viertel") == "viertel"


def _zelle(cfg) -> bool:
    return bool(cfg["field"].get("einheitszelle", False))


def _safe(name):
    keep = ("abcdefghijklmnopqrstuvwxyz"
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")
    return "".join(c if c in keep else "_" for c in str(name))


def schichtstapel(cfg):
    z = cfg["domain"]["z_basis_m"]
    out = []
    for L in reversed(list(cfg["layers"])):
        out.append({**L, "z_low": z, "z_high": z + float(L["maechtigkeit_m"])})
        z += float(L["maechtigkeit_m"])
    return out, z


def sondenpositionen(cfg) -> list[tuple[float, float]]:
    if _zelle(cfg):
        return [(0.0, 0.0)]
    nx, ny = int(cfg["field"]["n_x_full"]), int(cfg["field"]["n_y_full"])
    s = float(cfg["field"]["abstand_m"])
    xs = [-(nx - 1) * s / 2 + i * s for i in range(nx)]
    ys = [-(ny - 1) * s / 2 + i * s for i in range(ny)]
    if not _viertel(cfg):
        return [(x, y) for x in xs for y in ys]
    if nx % 2 or ny % 2:
        raise ValueError(
            f"Viertelmodell verlangt gerade Sondenzahlen je Richtung "
            f"(n_x_full={nx}, n_y_full={ny}); sonst läge eine Sonde auf einer "
            f"Symmetrieebene und müsste halbiert werden. Entweder Sondenzahl "
            f"ändern oder field.symmetrie = 'voll' setzen.")
    return [(x, y) for x in xs if x > 1e-9 for y in ys if y > 1e-9]


def temperaturleitfaehigkeit(cfg) -> float:
    lam = rc = ges = 0.0
    for L in cfg["layers"]:
        e = effektive_stoffwerte(L, cfg["fluid"])
        d = L["maechtigkeit_m"]
        lam += e["lambda"] * d; rc += e["rho_c"] * d; ges += d
    return (lam / ges) / (rc / ges)


def eindringtiefe_m(cfg) -> float:
    t = cfg["zeit"]["jahre"] * YEAR_D * DAY
    return float(np.sqrt(4.0 * temperaturleitfaehigkeit(cfg) * t))


def feldausdehnung(cfg):
    s = float(cfg["field"]["abstand_m"])
    return ((cfg["field"]["n_x_full"] - 1) * s / 2.0,
            (cfg["field"]["n_y_full"] - 1) * s / 2.0)


def gebiet_halbmasse(cfg):
    d = cfg["domain"]
    if not d.get("automatisch", True):
        return float(d["size_x_m"]), float(d["size_y_m"])
    fx, fy = feldausdehnung(cfg)
    rand = max(d.get("rand_min_m", 40.0),
               d.get("rand_faktor", 2.0) * eindringtiefe_m(cfg))
    return fx + rand, fy + rand


def gebiet(cfg):
    if _zelle(cfg):
        h = cfg["field"]["abstand_m"] / 2.0
        return -h, -h, h, h
    LX, LY = gebiet_halbmasse(cfg)
    return (0.0, 0.0, LX, LY) if _viertel(cfg) else (-LX, -LY, LX, LY)


def vertikale_intervalle(cfg):
    """Extrusionsintervalle: an jeder Schichtgrenze und an Sondenkopf/-fuß."""
    layers, z_top = schichtstapel(cfg)
    z_basis = cfg["domain"]["z_basis_m"]
    bh = cfg["borehole"]
    z_bh_top, z_bh_bot = z_top - bh["tiefe_kopf_m"], z_top - bh["tiefe_fuss_m"]
    if z_bh_bot < z_basis - 1e-9 or z_bh_top > z_top + 1e-9:
        raise ValueError("Sondentiefe liegt außerhalb der Schichtdomäne.")
    cuts = sorted({z_basis, z_top, z_bh_top, z_bh_bot}
                  | {L["z_low"] for L in layers} | {L["z_high"] for L in layers})
    mh = netzparameter(cfg)
    iv = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        if b - a < 1e-9:
            continue
        zc = 0.5 * (a + b)
        li = next(i for i, L in enumerate(layers)
                  if L["z_low"] - 1e-9 <= zc <= L["z_high"] + 1e-9)
        in_bh = (z_bh_bot - 1e-9) <= zc <= (z_bh_top + 1e-9)
        dz = mh["dz_sondenbereich_m"] if in_bh else mh["dz_ausserhalb_m"]
        iv.append({"z_low": a, "z_high": b, "layer": li, "in_bh": in_bh,
                   "n": max(1, int(np.ceil((b - a) / dz - 1e-9)))})
    return iv, layers, z_top, z_bh_bot, z_bh_top


# ======================================================================
#  Netz
# ======================================================================
def _grundriss(cfg, mh):
    """2D-Geometrie samt Größenfeldern. gmsh muss initialisiert sein."""
    import gmsh
    geo = gmsh.model.geo
    x_lo, y_lo, x_hi, y_hi = gebiet(cfg)
    z0 = cfg["domain"]["z_basis_m"]
    dx, dy = cfg["borehole"]["box_dx_m"], cfg["borehole"]["box_dy_m"]

    def rechteck(a, b, c, d):
        p = [geo.addPoint(a, b, z0), geo.addPoint(c, b, z0),
             geo.addPoint(c, d, z0), geo.addPoint(a, d, z0)]
        ln = [geo.addLine(p[i], p[(i + 1) % 4]) for i in range(4)]
        return geo.addCurveLoop(ln), ln

    cl_out, _ = rechteck(x_lo, y_lo, x_hi, y_hi)
    loops, kurven = [], []
    for x, y in sondenpositionen(cfg):
        cl, ln = rechteck(x - dx / 2, y - dy / 2, x + dx / 2, y + dy / 2)
        loops.append(cl); kurven += ln
    s_bg = geo.addPlaneSurface([cl_out] + loops)
    s_box = [geo.addPlaneSurface([c]) for c in loops]
    geo.synchronize()

    f_d = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_d, "CurvesList", sorted(set(kurven)))
    gmsh.model.mesh.field.setNumber(f_d, "Sampling", 8)
    f_1 = gmsh.model.mesh.field.add("Threshold")
    for k, v in [("InField", f_d), ("SizeMin", mh["groesse_an_sonde_m"]),
                 ("SizeMax", mh["groesse_im_feld_m"]), ("DistMin", 0.0),
                 ("DistMax", mh["ring_radius_m"])]:
        gmsh.model.mesh.field.setNumber(f_1, k, v)
    f_2 = gmsh.model.mesh.field.add("Threshold")
    for k, v in [("InField", f_d), ("SizeMin", mh["groesse_im_feld_m"]),
                 ("SizeMax", mh["groesse_fern_m"]),
                 ("DistMin", mh["feld_radius_m"]),
                 ("DistMax", mh["fern_radius_m"])]:
        gmsh.model.mesh.field.setNumber(f_2, k, v)
    f_m = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(f_m, "FieldsList", [f_1, f_2])
    gmsh.model.mesh.field.setAsBackgroundMesh(f_m)
    # Punktgrößen AUS: sonst schlägt die Größe an den Boxecken weit ins Feld
    # durch und groesse_im_feld_m wäre wirkungslos.
    for opt, val in [("Mesh.MeshSizeExtendFromBoundary", 0),
                     ("Mesh.MeshSizeFromPoints", 0),
                     ("Mesh.MeshSizeFromCurvature", 0),
                     ("Mesh.MshFileVersion", 2.2)]:
        gmsh.option.setNumber(opt, val)
    return [s_bg] + s_box, kurven


def knoten_2d(cfg, mh) -> int:
    """Knoten des 2D-Grundrisses — wirklich vernetzt, nicht geschätzt.

    Eine analytische Formel unterschätzt hier grob, weil der Übergang zwischen
    Feld und Fernfeld flächenmäßig dominiert. Die 2D-Vernetzung kostet
    Sekunden und ist exakt.
    """
    import gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("abzaehlen")
        _grundriss(cfg, mh)
        gmsh.model.mesh.generate(2)
        return len(gmsh.model.mesh.getNodes()[0])
    finally:
        gmsh.finalize()


_NETZ_CACHE: dict = {}


def netzparameter(cfg) -> dict:
    """Effektive Netzlängen — automatisch aus dem Feld oder von Hand."""
    nz = cfg["netz"]
    if not nz.get("automatisch", True):
        return {k: nz[k] for k in
                ("groesse_an_sonde_m", "ring_radius_m", "groesse_im_feld_m",
                 "feld_radius_m", "groesse_fern_m", "fern_radius_m",
                 "dz_sondenbereich_m", "dz_ausserhalb_m")} | {"vergroebert": 1.0}

    schluessel = (cfg["field"]["n_x_full"], cfg["field"]["n_y_full"],
                  cfg["field"]["abstand_m"], cfg["field"]["symmetrie"],
                  _zelle(cfg), cfg["zeit"]["jahre"], nz["ziel_knoten"])
    if schluessel in _NETZ_CACHE:
        return _NETZ_CACHE[schluessel]

    s = float(cfg["field"]["abstand_m"])
    b = float(cfg["borehole"]["box_dx_m"])
    L_bh = cfg["borehole"]["tiefe_fuss_m"] - cfg["borehole"]["tiefe_kopf_m"]
    dz_bh = L_bh / float(nz.get("dz_je_sonde", 40))
    n_lagen = int(np.ceil(L_bh / dz_bh)) + len(cfg["layers"]) + 1

    def satz(f):
        return {"groesse_an_sonde_m": b / float(nz["elemente_je_box"]) * f,
                "ring_radius_m": s / 3.0,
                "groesse_im_feld_m": s / 5.0 * f,
                "feld_radius_m": 1.5 * s,
                "groesse_fern_m": 3.0 * s * f,
                "fern_radius_m": 10.0 * s,
                "dz_sondenbereich_m": dz_bh,
                "dz_ausserhalb_m": 3.0 * dz_bh}

    ziel = float(nz["ziel_knoten"])
    f_max = float(nz["elemente_je_box"])      # nie gröber als die Boxbreite
    f = 1.0
    mh = satz(f)
    knoten = knoten_2d(cfg, mh) * n_lagen
    for _ in range(3):
        if knoten <= ziel or f >= f_max - 1e-9:
            break
        f = min(f * float(np.sqrt(knoten / ziel)), f_max)
        mh = satz(f)
        knoten = knoten_2d(cfg, mh) * n_lagen
    mh |= {"vergroebert": f, "knoten": knoten, "knotenlagen": n_lagen}
    _NETZ_CACHE[schluessel] = mh
    return mh


def netz_bauen(cfg: dict, out_dir: Path) -> Path:
    """2D-Grundriss Intervall für Intervall in z extrudieren (Prismen)."""
    import gmsh
    msh = out_dir / f"{cfg['ausgabe']['prefix']}.msh"
    x_lo, y_lo, x_hi, y_hi = gebiet(cfg)
    iv, layers, z_top, z_bh_bot, z_bh_top = vertikale_intervalle(cfg)
    dx, dy = cfg["borehole"]["box_dx_m"], cfg["borehole"]["box_dy_m"]
    for x, y in sondenpositionen(cfg):
        if not (x_lo < x - dx / 2 and x + dx / 2 < x_hi
                and y_lo < y - dy / 2 and y + dy / 2 < y_hi):
            raise ValueError(f"Sondenbox bei ({x:.1f}, {y:.1f}) passt nicht ins "
                             f"Gebiet — domain vergrößern.")
    mh = netzparameter(cfg)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("btes")
    geo = gmsh.model.geo
    basis, kurven = _grundriss(cfg, mh)

    # Eine Extrusion über alles auf einmal ergäbe nur EIN Volumen je
    # Grundfläche; für eigene Materialien je Schicht braucht es getrennte.
    vol_layer = {i: [] for i in range(len(layers))}
    vol_bh, cur = [], list(basis)
    for it in iv:
        ext = geo.extrude([(2, t) for t in cur], 0, 0,
                          it["z_high"] - it["z_low"],
                          numElements=[it["n"]], recombine=True)
        geo.synchronize()
        tops, vols, prev = [], [], None
        for d, t in ext:
            if d == 2:
                prev = t
            elif d == 3:
                vols.append(t); tops.append(prev)
        vol_layer[it["layer"]].append(vols[0])
        (vol_bh if it["in_bh"] else vol_layer[it["layer"]]).extend(vols[1:])
        cur = tops

    surf_bot, surf_top = list(basis), list(cur)
    surf_aussen, surf_sym = [], []

    def flach(lo, hi, v):
        return abs(lo - v) < 1e-6 and abs(hi - v) < 1e-6

    for _, tag in gmsh.model.getEntities(2):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
        if abs(zmax - zmin) < 1e-6:
            continue
        sym = (_viertel(cfg) and not _zelle(cfg)
               and (flach(xmin, xmax, 0.0) or flach(ymin, ymax, 0.0)))
        aus = (flach(xmin, xmax, x_hi) or flach(ymin, ymax, y_hi)
               or flach(xmin, xmax, x_lo) or flach(ymin, ymax, y_lo))
        (surf_sym if sym else surf_aussen if aus else []).append(tag)

    pg = 1
    for i, L in enumerate(layers):
        gmsh.model.addPhysicalGroup(3, vol_layer[i], tag=pg, name=L["name"]); pg += 1
    gmsh.model.addPhysicalGroup(3, vol_bh, tag=pg, name="sonden")
    gmsh.model.addPhysicalGroup(2, surf_top, tag=200, name="top")
    gmsh.model.addPhysicalGroup(2, surf_bot, tag=201, name="bottom")
    if surf_aussen:
        gmsh.model.addPhysicalGroup(2, surf_aussen, tag=202, name="lateral")
    if surf_sym:
        # bewusst OHNE Randbedingung: eine Fläche ohne Randbedingung ist in
        # OGS ein Nullfluss-Rand, und das ist genau die Symmetriebedingung.
        gmsh.model.addPhysicalGroup(2, surf_sym, tag=203, name="symmetrie")

    gmsh.model.mesh.generate(3)
    gmsh.write(str(msh))
    n = len(gmsh.model.mesh.getNodes()[0])
    gmsh.finalize()
    print(f"      {len(sondenpositionen(cfg))} Sonden, {len(iv)} Intervalle, "
          f"{sum(i['n'] for i in iv)} Elementlagen, {n:,} Knoten")
    return msh


def netz_dateien(cfg) -> dict:
    p = cfg["ausgabe"]["prefix"]
    return {"domain": f"{p}_domain.vtu", "top": f"{p}_physical_group_top.vtu",
            "bottom": f"{p}_physical_group_bottom.vtu",
            "sonden": f"{p}_physical_group_sonden.vtu"}


def netz_wandeln(cfg, msh, out_dir) -> dict:
    import ogstools as ot
    meshes = ot.Meshes.from_gmsh(filename=str(msh), dim=3, reindex=True, log=False)
    pref = cfg["ausgabe"]["prefix"]
    for name, mesh in meshes.items():
        fn = (f"{pref}_domain.vtu" if name == "domain"
              else f"{pref}_physical_group_{_safe(name)}.vtu")
        mesh.save(str(Path(out_dir) / fn), binary=True)
    return netz_dateien(cfg)


# ======================================================================
#  Lastkurve und Projektdatei
# ======================================================================
def lastkurve(cfg) -> dict:
    P = leistung_je_sonde_W(cfg)
    P_ref = float(np.abs(P).max())
    if P_ref == 0:
        raise ValueError("Lastprofil ist überall null.")
    rampe = max(60.0, cfg["last"]["rampe_tage"] * DAY)
    monat = MONAT_D * DAY
    t, v, now = [0.0], [0.0], 0.0
    for _ in range(int(cfg["zeit"]["jahre"])):
        for p in P:
            q = float(p) / P_ref
            now += rampe; t.append(now); v.append(q)
            halt = max(0.0, monat - rampe)
            if halt > 0:
                now += halt; t.append(now); v.append(q)
    now += rampe; t.append(now); v.append(0.0)
    return {"t_end": now, "kurve": (np.array(t), np.array(v)), "P_ref_W": P_ref}


def _se(parent, tag, text=None, **attrs):
    el = ET.SubElement(parent, tag, **{k: str(v) for k, v in attrs.items()})
    if text is not None:
        el.text = str(text)
    return el


def _prop(parent, name, value):
    p = _se(parent, "property")
    _se(p, "name", name); _se(p, "type", "Constant"); _se(p, "value", value)


def _indent(el, lvl=0):
    pad = "\n" + lvl * "    "
    if len(el):
        if not (el.text and el.text.strip()):
            el.text = pad + "    "
        for c in el:
            _indent(c, lvl + 1)
        if not (el[-1].tail and el[-1].tail.strip()):
            el[-1].tail = pad
    if lvl and not (el.tail and el.tail.strip()):
        el.tail = pad


def prj_bauen(cfg, out_dir: Path, nf: dict, kv: dict) -> Path:
    pref, sol, fluid = cfg["ausgabe"]["prefix"], cfg["loeser"], cfg["fluid"]
    bh = cfg["borehole"]
    V_box = bh["box_dx_m"] * bh["box_dy_m"] * (bh["tiefe_fuss_m"]
                                               - bh["tiefe_kopf_m"])
    q_v, T0 = kv["P_ref_W"] / V_box, cfg["initial"]["T_C"] + 273.15

    r = ET.Element("OpenGeoSysProject")
    m = _se(r, "meshes")
    for k in ("domain", "top", "bottom", "sonden"):
        _se(m, "mesh", nf[k])

    p = _se(_se(r, "processes"), "process")
    _se(p, "name", "HeatConduction"); _se(p, "type", "HEAT_CONDUCTION")
    _se(p, "integration_order", 2)
    _se(_se(p, "process_variables"), "process_variable", "T")

    media = _se(r, "media")
    lb = list(reversed(cfg["layers"]))
    for i, L in enumerate(lb):
        e = effektive_stoffwerte(L, fluid)
        med = _se(media, "medium", id=i); _se(med, "phases")
        pr = _se(med, "properties")
        _prop(pr, "thermal_conductivity", e["lambda"])
        _prop(pr, "density", e["dichte"])
        _prop(pr, "specific_heat_capacity", e["cp"])
    # Sondenmaterial: mittlere Schicht (didaktische Voreinstellung; für eine
    # realistische Rechnung gehört hier ein Verfüllmaterial hinein)
    e = effektive_stoffwerte(lb[len(lb) // 2], fluid)
    med = _se(media, "medium", id=len(lb)); _se(med, "phases")
    pr = _se(med, "properties")
    _prop(pr, "thermal_conductivity", e["lambda"])
    _prop(pr, "density", e["dichte"])
    _prop(pr, "specific_heat_capacity", e["cp"])

    tl = _se(r, "time_loop")
    pr_ = _se(_se(tl, "processes"), "process", ref="HeatConduction")
    _se(pr_, "nonlinear_solver", "picard")
    cc = _se(pr_, "convergence_criterion")
    _se(cc, "type", "DeltaX"); _se(cc, "norm_type", "NORM2")
    _se(cc, "reltol", sol["reltol_T"])
    _se(_se(pr_, "time_discretization"), "type", "BackwardEuler")
    ts = _se(pr_, "time_stepping")
    _se(ts, "type", "FixedTimeStepping")
    _se(ts, "t_initial", 0.0); _se(ts, "t_end", kv["t_end"])
    dt = MONAT_D * DAY / cfg["zeit"]["schritte_je_monat"]
    n_st = int(np.ceil(kv["t_end"] / dt))
    pa = _se(_se(ts, "timesteps"), "pair")
    _se(pa, "repeat", n_st); _se(pa, "delta_t", dt)

    o = _se(tl, "output")
    _se(o, "type", "VTK"); _se(o, "prefix", pref)
    pa = _se(_se(o, "timesteps"), "pair")
    _se(pa, "repeat", n_st)
    _se(pa, "each_steps", cfg["zeit"]["ausgabe_je_n_schritte"])
    _se(o, "output_iteration_results", "false")
    _se(_se(o, "variables"), "variable", "T")

    par = _se(r, "parameters")
    g = cfg["initial"]["geothermischer_gradient_K_m"]
    e_ = _se(par, "parameter"); _se(e_, "name", "T0")
    if abs(g) > 1e-12:
        z_tot = sum(L["maechtigkeit_m"] for L in cfg["layers"])
        _se(e_, "type", "Function")
        _se(e_, "expression", f"{T0} + ({g:.6g})*({z_tot} - z)")
    else:
        _se(e_, "type", "Constant"); _se(e_, "value", T0)
    e_ = _se(par, "parameter")
    _se(e_, "name", "q_amp"); _se(e_, "type", "Constant"); _se(e_, "value", q_v)
    e_ = _se(par, "parameter")
    _se(e_, "name", "q_sonde"); _se(e_, "type", "CurveScaled")
    _se(e_, "curve", "last"); _se(e_, "parameter", "q_amp")

    c = _se(_se(r, "curves"), "curve")
    tt, vv = kv["kurve"]
    _se(c, "name", "last")
    _se(c, "coords", " ".join(f"{x:.6e}" for x in tt))
    _se(c, "values", " ".join(f"{x:.6e}" for x in vv))

    pv = _se(_se(r, "process_variables"), "process_variable")
    _se(pv, "name", "T"); _se(pv, "components", 1); _se(pv, "order", 1)
    _se(pv, "initial_condition", "T0")
    bcs = _se(pv, "boundary_conditions")
    for f in ("top", "bottom"):
        b = _se(bcs, "boundary_condition")
        _se(b, "mesh", Path(nf[f]).stem)
        _se(b, "type", "Dirichlet"); _se(b, "parameter", "T0")
    # Seitenflächen bleiben ohne Randbedingung: auf Symmetrieebenen exakt die
    # Spiegelbedingung, außen adiabat.
    st = _se(_se(pv, "source_terms"), "source_term")
    _se(st, "mesh", Path(nf["sonden"]).stem)
    _se(st, "type", "Volumetric"); _se(st, "parameter", "q_sonde")

    nl = _se(_se(r, "nonlinear_solvers"), "nonlinear_solver")
    _se(nl, "name", "picard"); _se(nl, "type", "Picard")
    _se(nl, "max_iter", sol["nichtlinear_max_iter"])
    _se(nl, "linear_solver", "ls")
    ls = _se(_se(r, "linear_solvers"), "linear_solver")
    _se(ls, "name", "ls")
    ei = _se(ls, "eigen")
    _se(ei, "solver_type", sol["typ"])
    _se(ei, "precon_type", sol["vorkonditionierer"])
    _se(ei, "max_iteration_step", sol["max_iter"])
    _se(ei, "error_tolerance", sol["toleranz"])
    _se(ei, "scaling", "true")

    _indent(r)
    prj = out_dir / f"{pref}.prj"
    ET.ElementTree(r).write(prj, encoding="ISO-8859-1", xml_declaration=True)
    return prj


def ogs_starten(prj: Path) -> int:
    exe = shutil.which("ogs") or shutil.which("ogs.exe")
    if not exe:
        print("ogs nicht im PATH", file=sys.stderr)
        return 1
    print(">>", exe, prj, "-o", prj.parent)
    return subprocess.call([exe, str(prj), "-o", str(prj.parent)])


# ======================================================================
#  Bericht
# ======================================================================
def bericht(cfg) -> None:
    lp = lastprofil(cfg)
    n_full, n_mod = _n_full(cfg), len(sondenpositionen(cfg))
    P_bh = leistung_je_sonde_W(cfg)
    L_bh = cfg["borehole"]["tiefe_fuss_m"] - cfg["borehole"]["tiefe_kopf_m"]
    modus = ("Einheitszelle" if _zelle(cfg)
             else "Viertelmodell" if _viertel(cfg) else "Vollmodell")

    print("\n--- Last ---")
    if "solar_MWh" in lp:
        print(f"  Solarertrag        {lp['solar_MWh']:>10,.0f} MWh/a "
              f"({cfg['last']['kollektorflaeche_m2']:,.0f} m² Kollektor)")
        print(f"  Wärmebedarf        {lp['bedarf_MWh']:>10,.0f} MWh/a "
              f"(Reduktionsgrad "
              f"{cfg['last'].get('reduktionsgrad_prozent', 0.0):.1f} %)")
    print(f"  Beladung           {lp['beladung_MWh']:>10,.0f} MWh/a")
    print(f"  Entnahme           {lp['entnahme_MWh']:>10,.0f} MWh/a   = "
          f"{100*lp['entnahme_MWh']/max(lp['beladung_MWh'],1e-9):.0f} % der Beladung")
    print(f"  Jahresbilanz       {lp['bilanz_MWh']:>+10,.0f} MWh/a", end="")
    print("   <-- Speicher wird netto entleert, das Feld kühlt aus"
          if lp["bilanz_MWh"] < -1 else
          "   <-- Überschuss, der nicht abgenommen wird"
          if lp["bilanz_MWh"] > 1 else "   <-- ausgeglichen")
    if "reduktion_bilanz_null" in lp:
        print(f"  Reduktionsgrad für ausgeglichene Bilanz: "
              f"{lp['reduktion_bilanz_null']:.1f} %")

    print(f"\n--- Feld ({modus}) ---")
    print(f"  Sonden gesamt {n_full}, davon {n_mod} gerechnet")
    print(f"  Spitzenlast je Sonde {np.abs(P_bh).max():,.0f} W = "
          f"{np.abs(P_bh).max()/L_bh:.1f} W/m "
          f"(Literatur Erdwärmesonden 20–70 W/m)")
    print(f"\n  {'Monat':<6}{'Feld [kW]':>12}{'je Sonde [W]':>14}{'W/m':>9}")
    for i, mn in enumerate(MONATE):
        print(f"  {mn:<6}{lp['P_kW'][i]:>12,.1f}{P_bh[i]:>14,.0f}"
              f"{P_bh[i]/L_bh:>9.1f}")

    if not _zelle(cfg):
        fx, fy = feldausdehnung(cfg)
        LX, LY = gebiet_halbmasse(cfg)
        print("\n--- Gebiet ---")
        print(f"  Feld            {2*fx:>8,.0f} × {2*fy:<8,.0f} m")
        print(f"  Eindringtiefe   {eindringtiefe_m(cfg):>8,.0f} m   "
              f"(alpha = {temperaturleitfaehigkeit(cfg):.2e} m²/s, "
              f"{cfg['zeit']['jahre']} a)")
        print(f"  Gebiet          {2*LX:>8,.0f} × {2*LY:<8,.0f} m   "
              f"Rand {LX-fx:,.0f} m hinter der äußersten Sonde")

    mh = netzparameter(cfg)
    print("\n--- Netz ---")
    print(f"  an der Sondenwand {mh['groesse_an_sonde_m']:>6.2f} m   "
          f"im Feld {mh['groesse_im_feld_m']:>6.2f} m   "
          f"fern {mh['groesse_fern_m']:>6.1f} m")
    print(f"  Elementhöhe       {mh['dz_sondenbereich_m']:>6.2f} m im "
          f"Sondenbereich, {mh['dz_ausserhalb_m']:.1f} m darüber/darunter")
    if "knoten" in mh:
        print(f"  Knoten            {mh['knoten']:>10,.0f}   = Freiheitsgrade "
              f"(nur T), Budget {cfg['netz']['ziel_knoten']:,}")
        if mh["vergroebert"] > 1.001:
            print(f"  ACHTUNG: Ebene um Faktor {mh['vergroebert']:.2f} "
                  f"vergröbert, damit das Budget hält.")
            print("           Das Ergebnis liegt dadurch auf der optimistischen "
                  "Seite — ein zu grobes Netz macht die Sonde zu warm.")
            print("           Gegenmittel: netz.ziel_knoten erhöhen, weniger "
                  "Sonden, oder field.symmetrie = 'viertel'.")


# ======================================================================
#  Darstellung
# ======================================================================
SURF, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"


def de(x, nd=1):
    """Deutsche Zahlschreibweise."""
    s = f"{x:,.{nd}f}"
    return s.replace(",", "\x00").replace(".", ",").replace("\x00", ".")


def _plt():
    """matplotlib mit dem Stil dieser Auswertung."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": SURF, "axes.facecolor": SURF,
        "savefig.facecolor": SURF, "text.color": INK,
        "axes.labelcolor": INK2, "axes.edgecolor": AXIS,
        "xtick.color": MUTED, "ytick.color": MUTED, "grid.color": GRID,
        "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10,
        "axes.titleweight": "bold", "axes.titlelocation": "left",
        "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 1.6, "legend.frameon": False})
    return plt


def _fmt():
    from matplotlib.ticker import FuncFormatter
    return FuncFormatter(lambda v, _: de(v, 0))


def _figdir(cfg) -> Path:
    d = HERE / cfg["ausgabe"]["figures"]
    d.mkdir(exist_ok=True)
    return d


# ======================================================================
#  Auswertung des Laufs
# ======================================================================
def zeitreihe(cfg, out: Path):
    """Zeitreihe aus den Ergebnisdateien. Sondenlage und Last kommen aus dem
    Lauf selbst, nicht aus CONFIG — sonst passt beides nicht zusammen, wenn
    die Einstellungen nach dem Lauf geändert wurden."""
    import pyvista as pv
    pref = cfg["ausgabe"]["prefix"]
    files = sorted(glob.glob(str(out / f"{pref}_ts_*.vtu")),
                   key=lambda f: float(re.search(r"_t_([\d.]+)\.vtu$", f).group(1)))
    if not files:
        raise SystemExit(f"Keine Ergebnisdateien in {out}.")
    dom = pv.read(str(out / f"{pref}_domain.vtu"))
    mid = np.asarray(dom.cell_data["MaterialIDs"])
    vol = np.abs(np.asarray(dom.compute_cell_sizes(
        length=False, area=False, volume=True).cell_data["Volume"]))
    cc = np.asarray(dom.cell_centers().points)
    pts = np.asarray(dom.points)
    n_lay = len(cfg["layers"])
    T0K = cfg["initial"]["T_C"] + 273.15

    rho_c = np.zeros(dom.n_cells)
    lb = list(reversed(cfg["layers"]))
    for i, L in enumerate(lb):
        rho_c[mid == i] = effektive_stoffwerte(L, cfg["fluid"])["rho_c"]
    rho_c[mid >= n_lay] = effektive_stoffwerte(
        lb[len(lb) // 2], cfg["fluid"])["rho_c"]

    dxy = cfg["borehole"]["box_dx_m"]
    is_bh = mid >= n_lay
    if not is_bh.any():
        raise SystemExit(f"Keine Sondenzellen in {out} (MaterialID >= {n_lay}).")

    def cluster(v, tol):
        v = np.sort(np.unique(np.round(v, 3)))
        gr, cur = [], [v[0]]
        for x in v[1:]:
            if x - cur[-1] <= tol:
                cur.append(x)
            else:
                gr.append(float(np.mean(cur))); cur = [x]
        gr.append(float(np.mean(cur)))
        return gr

    masken, positionen = [], []
    for x in cluster(cc[is_bh, 0], dxy):
        for y in cluster(cc[is_bh, 1], dxy):
            m = (is_bh & (np.abs(cc[:, 0] - x) < dxy)
                 & (np.abs(cc[:, 1] - y) < dxy))
            if m.sum():
                masken.append(m); positionen.append((x, y))
    gew = [vol[m] / vol[m].sum() for m in masken]

    x_lo, y_lo, x_hi, y_hi = gebiet(cfg)
    tol = 1e-6
    rand = ((np.abs(pts[:, 0] - x_hi) < tol) | (np.abs(pts[:, 1] - y_hi) < tol)
            | (np.abs(pts[:, 0] - x_lo) < tol) | (np.abs(pts[:, 1] - y_lo) < tol))
    if _viertel(cfg) and not _zelle(cfg):
        rand &= ~((np.abs(pts[:, 0]) < tol) | (np.abs(pts[:, 1]) < tol))

    wurzel = ET.parse(str(out / f"{pref}.prj")).getroot()
    # Der Name der Quellstärke hat sich über die Fassungen geändert; hier
    # werden beide akzeptiert, damit auch ältere Läufe auswertbar bleiben.
    q_amp = next((float(p.findtext("value")) for p in wurzel.iter("parameter")
                  if p.findtext("name") in ("q_amp", "q_v_amp")), None)
    if q_amp is None:
        namen = [p.findtext("name") for p in wurzel.iter("parameter")]
        raise SystemExit(
            f"In {pref}.prj findet sich keine Quellstärke (erwartet 'q_amp' "
            f"oder 'q_v_amp'). Vorhandene Parameter: {namen}")
    kv = wurzel.find("curves/curve")
    k_t = np.array([float(x) for x in kv.findtext("coords").split()])
    k_v = np.array([float(x) for x in kv.findtext("values").split()])
    A_box = cfg["borehole"]["box_dx_m"] * cfg["borehole"]["box_dy_m"]
    korr = fluid_korrektur(cfg)

    print(f"{len(files)} Ausgabezeitpunkte, {len(masken)} Sonden, "
          f"{rand.sum():,} Randknoten")
    print(f"Last aus {pref}.prj: Förderspitze {k_v.min()*q_amp*A_box:.1f} W/m, "
          f"Fluidkorrektur {k_v.min()*q_amp*A_box*korr:+.2f} K")

    zeilen, je_sonde = [], []
    for i, f in enumerate(files):
        td = float(re.search(r"_t_([\d.]+)\.vtu$", f).group(1)) / DAY
        m = pv.read(f)
        Tp = np.asarray(m.point_data["T"])
        Tc = np.asarray(m.point_data_to_cell_data().cell_data["T"])
        bm = np.array([float((Tc[k] * w).sum()) for k, w in zip(masken, gew)])
        q = q_amp * float(np.interp(td * DAY, k_t, k_v)) * A_box
        je_sonde.append(bm - 273.15)
        zeilen.append(dict(
            t_a=td / YEAR_D,
            T_wand_min_C=bm.min() - 273.15, T_wand_mittel_C=bm.mean() - 273.15,
            T_wand_max_C=bm.max() - 273.15,
            T_fluid_min_C=bm.min() - 273.15 + q * korr, q_W_m=q,
            E_GJ=float((rho_c * (Tc - T0K) * vol).sum()) / 1e9,
            dT_rand_K=float(np.abs(Tp[rand] - T0K).max()) if rand.any() else 0.0))
        if (i + 1) % 60 == 0 or i == len(files) - 1:
            print(f"  {i+1}/{len(files)}  t = {td/YEAR_D:5.2f} a  "
                  f"Fluid {zeilen[-1]['T_fluid_min_C']:+7.2f} °C")
    return zeilen, positionen, np.array(je_sonde), korr


def auswerten(cfg, out: Path) -> None:
    """Vereisungsprüfung: Tabellen, Abbildung und Urteil im Klartext.

    Geschrieben werden drei Tabellen:
      summary.csv            je Zeitschritt die Kennwerte über alle Sonden
      sonden_temperaturen.csv  je Zeitschritt JEDE Sonde einzeln
      sonden_positionen.csv    Lage und Kennung der Sonden
    """
    z, positionen, je_sonde, korr = zeitreihe(cfg, out)
    with open(HERE / "summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(z[0]))
        w.writeheader(); w.writerows(z)
    print(f"\n-> {HERE / 'summary.csv'}")

    # Lage der Sonden, aus dem Netz gelesen
    with open(HERE / "sonden_positionen.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sonde", "x_m", "y_m"])
        for i, (x, y) in enumerate(positionen):
            w.writerow([f"S{i:03d}", f"{x:.3f}", f"{y:.3f}"])
    print(f"-> {HERE / 'sonden_positionen.csv'}  ({len(positionen)} Sonden)")

    # Jede Sonde einzeln über die Zeit. Spalten S000, S001, ... in der
    # Reihenfolge der Positionstabelle; Werte sind die volumengewichtete
    # Temperatur der Sondenbox in °C. Die Fluidtemperatur folgt daraus als
    #     T_Fluid = T_Sonde + q_W_m * korr
    # mit q_W_m aus derselben Zeile und korr = fluid_korrektur(CONFIG).
    with open(HERE / "sonden_temperaturen.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["t_a", "q_W_m", "fluid_korrektur_K_pro_W_m"]
                   + [f"S{i:03d}" for i in range(len(positionen))])
        for zeile, temps in zip(z, je_sonde):
            w.writerow([f"{zeile['t_a']:.6f}", f"{zeile['q_W_m']:.4f}",
                        f"{korr:.6f}"] + [f"{v:.4f}" for v in temps])
    print(f"-> {HERE / 'sonden_temperaturen.csv'}  "
          f"({len(z)} Zeitpunkte × {len(positionen)} Sonden)")

    t = np.array([r["t_a"] for r in z])
    fl = np.array([r["T_fluid_min_C"] for r in z])
    wand = np.array([r["T_wand_min_C"] for r in z])
    E = np.array([r["E_GJ"] for r in z])
    dR = np.array([r["dT_rand_K"] for r in z])
    jahre = np.arange(1, max(1, int(np.floor(t.max()))) + 1)
    jmin = np.array([fl[(t > j - 1) & (t <= j)].min()
                     if ((t > j - 1) & (t <= j)).any() else fl.min()
                     for j in jahre])

    plt = _plt(); F = _fmt()
    fig, ax = plt.subplots(1, 3, figsize=(14.0, 4.5),
                           gridspec_kw={"width_ratios": [1.25, 1.0, 1.0]})
    # Fluidtemperatur jeder einzelnen Sonde; als Band die Spreizung über das
    # Feld, als kräftige Linie die kälteste Sonde — sie entscheidet.
    q_t = np.array([r["q_W_m"] for r in z])
    fl_alle = je_sonde + (q_t * korr)[:, None]
    fl_kalt, fl_warm = fl_alle.min(axis=1), fl_alle.max(axis=1)

    n_s = je_sonde.shape[1]
    a = ax[0]
    a.axhline(cfg["initial"]["T_C"], color=MUTED, lw=0.8, ls=":")
    a.fill_between(t, fl_kalt, fl_warm, color=C1, alpha=0.16, lw=0)
    a.plot(t, wand, color=C1, lw=0.8, alpha=0.40)
    a.plot(t, fl_kalt, color=C1, lw=1.2)
    for (g, lab), col in zip(GRENZEN, (C4, C2)):
        a.axhline(g, color=col, lw=1.3, ls="--")
    j = t > t.max() - 1.0
    # Bei einer einzelnen Sonde sind kälteste und wärmste dieselbe; dann nur
    # Fluid und Wand beschriften, sonst stehen zwei Namen auf einer Linie.
    a.annotate("Fluid" if n_s == 1 else "kälteste Sonde\nFluid",
               (t.max(), fl_kalt[j].min()), color=C1,
               fontsize=8.5, fontweight="bold", xytext=(5, 0),
               textcoords="offset points", va="center")
    if n_s > 1:
        a.annotate("wärmste Sonde", (t.max(), fl_warm[j].max()), color=C1,
                   alpha=0.75, fontsize=8, xytext=(5, 0),
                   textcoords="offset points", va="center")
    a.annotate("Sondenwand", (t.max(), wand[j].max()), color=C1, alpha=0.6,
               fontsize=8, xytext=(5, 0), textcoords="offset points",
               va="center")
    # Die Grenzlinien liegen nur 4 K auseinander und laufen quer durch das
    # Band — beschriftet werden sie deshalb in Panel 2, wo Platz ist.
    # Rechts nur so viel Rand, wie die Beschriftung braucht.
    a.set_xlim(0, t.max() * 1.17)
    a.set_xlabel("Betriebsjahr"); a.set_ylabel("Temperatur [°C]")
    a.set_title("Fluidtemperatur der Einheitszelle — gestrichelt die Grenzwerte"
                if n_s == 1 else
                f"Fluidtemperatur aller {n_s} Sonden — gestrichelt die Grenzwerte")
    a.xaxis.set_major_formatter(F); a.yaxis.set_major_formatter(F)
    a.grid(axis="y", alpha=0.9); a.set_axisbelow(True)

    b = ax[1]
    for (g, lab), col in zip(GRENZEN, (C4, C2)):
        b.axhline(g, color=col, lw=1.2, ls="--")
        b.annotate(lab, (jahre[-1], g), color=col, fontsize=8.5,
                   fontweight="bold", ha="right", va="bottom",
                   xytext=(0, 4), textcoords="offset points")
    b.plot(jahre, jmin, color=C1, marker="o", ms=3.4, mfc=SURF, mew=1.2)
    # Kennzahl in die freie Ecke, nicht an die Kurve. Eine fallende Kurve
    # läuft von oben links nach unten rechts — frei ist dann unten links.
    # Bei fallender Kurve liegen ausserdem die Grenzlinien oben, wo die
    # Beschriftung sonst mit ihnen kollidiert.
    steigt = jmin[-1] > jmin[0]
    b.text(0.03, 0.94 if steigt else 0.06,
           f"Jahr {jahre[-1]}:  {de(jmin[-1], 2)} °C",
           transform=b.transAxes, color=C1, fontsize=9, fontweight="bold",
           va="top" if steigt else "bottom")
    b.set_xlabel("Betriebsjahr"); b.set_ylabel("Fluidminimum [°C]")
    b.set_title("Prüfung je Betriebsjahr")
    b.xaxis.set_major_formatter(F); b.yaxis.set_major_formatter(F)
    b.grid(axis="y", alpha=0.9); b.set_axisbelow(True)

    c = ax[2]
    c.plot(t, E, color=C3); c.axhline(0, color=AXIS, lw=1.0)
    c.set_xlabel("Betriebsjahr"); c.set_ylabel("Wärmeinhalt [GJ]")
    c.set_title(f"Wärmeinhalt — Rand max {de(dR.max(), 4)} K")
    c.xaxis.set_major_formatter(F); c.yaxis.set_major_formatter(F)
    c.grid(axis="y", alpha=0.9); c.set_axisbelow(True)

    friert = jmin.min() < 0.0
    fig.suptitle(f"Vereisungsprüfung — "
                 f"{'FRIERT EIN' if friert else 'friert nicht ein'}, tiefstes "
                 f"Fluidminimum {de(jmin.min(), 2)} °C", x=0.006, ha="left",
                 fontsize=11.5, fontweight="bold", color=C4 if friert else INK)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    ziel = _figdir(cfg) / "5_vereisung.png"
    fig.savefig(ziel, dpi=150); plt.close(fig)
    print(f"-> {ziel}")

    print("\n" + "=" * 62)
    print("VEREISUNGSPRÜFUNG")
    print("=" * 62)
    print(f"  simulierte Zeit        {t.max():>8.1f} a")
    print(f"  tiefstes Fluidminimum  {jmin.min():>8.2f} °C")
    print(f"  Sondenwand dazu        {wand.min():>8.2f} °C")
    for g, lab in GRENZEN:
        unter = jahre[jmin < g]
        print(f"  {lab:<20} " + (f"unterschritten ab Jahr {unter[0]}"
                                 if len(unter) else
                                 f"gehalten, Reserve {jmin.min()-g:+.2f} K"))
    if len(jmin) > 11:
        drift = jmin[-1] - jmin[-11]
        print(f"  Drift letzte 10 Jahre  {drift:>+8.2f} K", end="")
        print("   -> eingeschwungen" if abs(drift) < 0.3
              else "   -> läuft weiter weg")
    # Bei der Einheitszelle SIND die Seitenränder Symmetrieebenen des
    # unendlichen Feldes — dort gehört eine Störung hin und ist kein Fehler.
    if _zelle(cfg):
        print(f"  |dT| am Zellrand       {dR.max():>8.4f} K   -> "
              f"Symmetrieebene, Störung dort erwartet")
    else:
        print(f"  |dT| am Aussenrand     {dR.max():>8.4f} K", end="")
        print("   -> Gebiet gross genug" if dR.max() < 0.5
              else "   -> GEBIET ZU KLEIN, Rand faelscht mit")
    # Auslegung aus der Linearität — ohne weitere Rechnung
    print("\n  nötige Sondenzahl bei diesem Lastprofil:")
    for g, lab in GRENZEN:
        n = sondenzahl_fuer_grenze(cfg, float(jmin.min()), g)
        print(f"    {lab:<20} {n:>10,.0f}   (vorhanden {_n_full(cfg)})")


# ======================================================================
#  Abbildungen zur Konfiguration
# ======================================================================
def abbildungen(cfg) -> None:
    from matplotlib.patches import Rectangle
    from matplotlib.colors import LinearSegmentedColormap
    plt = _plt(); F = _fmt(); FIG = _figdir(cfg)
    SEQ = LinearSegmentedColormap.from_list(
        "b", ["#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#184f95", "#0d366b"])

    # --- 1: Feld und Symmetrie ---
    nx, ny = cfg["field"]["n_x_full"], cfg["field"]["n_y_full"]
    s = cfg["field"]["abstand_m"]
    xs = np.array([-(nx - 1) * s / 2 + i * s for i in range(nx)])
    ys = np.array([-(ny - 1) * s / 2 + i * s for i in range(ny)])
    X, Y = np.meshgrid(xs, ys)
    LX, LY = gebiet_halbmasse(cfg)
    fig, ax = plt.subplots(figsize=(9.6, 7.0))
    ax.add_patch(Rectangle((-LX, -LY), 2 * LX, 2 * LY, facecolor="#f4f3ef",
                           edgecolor=AXIS, lw=1.2, ls="--"))
    if _viertel(cfg) and not _zelle(cfg):
        ax.add_patch(Rectangle((0, 0), LX, LY, facecolor="#e8f0fb",
                               edgecolor=C1, lw=1.8))
        q = (X > 0) & (Y > 0)
        ax.plot(X[~q], Y[~q], "o", ms=3.2, color=MUTED, alpha=0.40)
        ax.plot(X[q], Y[q], "o", ms=3.4, color=C1)
        ax.axhline(0, color=C2, lw=1.9); ax.axvline(0, color=C2, lw=1.9)
        ax.annotate(f"gerechnet: {int(q.sum())} von {nx*ny} Sonden",
                    (LX * 0.55, LY * 0.66), color=C1, fontsize=9.5,
                    fontweight="bold", ha="center")
        ax.annotate("gespiegelt, nicht gerechnet", (-LX * 0.55, -LY * 0.66),
                    color=MUTED, fontsize=8.5, ha="center")
        ax.annotate("Symmetrieebenen x = 0 und y = 0\nkeine Randbedingung —\n"
                    "in OGS ist das exakt Nullfluss", (-LX * 0.55, LY * 0.62),
                    color=C2, fontsize=8.8, fontweight="bold", ha="center")
        titel = f"Sondenfeld und Viertelmodell — {nx} × {ny}, Raster {de(s,0)} m"
    else:
        ax.plot(X, Y, "o", ms=3.4, color=C1)
        ax.annotate(f"gerechnet: alle {nx*ny} Sonden", (0, LY * 0.72),
                    color=C1, fontsize=9.5, fontweight="bold", ha="center")
        titel = f"Sondenfeld, Vollmodell — {nx} × {ny}, Raster {de(s,0)} m"
    ax.annotate(f"Gebiet {de(2*LX,0)} × {de(2*LY,0)} m",
                (-LX * 0.96, -LY * 0.94), color=INK2, fontsize=8.5)
    ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    ax.xaxis.set_major_formatter(F); ax.yaxis.set_major_formatter(F)
    ax.set_title(titel, loc="left")
    for sp in ax.spines.values():
        sp.set_visible(False)
    fig.tight_layout(); fig.savefig(FIG / "1_feld.png", dpi=150); plt.close(fig)
    print("-> 1_feld.png")

    # --- 2: Untergrund ---
    layers = list(cfg["layers"])
    eff = [effektive_stoffwerte(L, cfg["fluid"]) for L in layers]
    lam = np.array([e["lambda"] for e in eff])
    z_top = sum(L["maechtigkeit_m"] for L in layers)
    norm = plt.Normalize(lam.min(), lam.max())
    fig, ax = plt.subplots(figsize=(11.0, 9.8))
    X0, X1, XT = 0.0, 0.20, 0.30
    rows = np.linspace(3.0, z_top - 3.0, len(layers))
    d = 0.0
    for i, (L, lm) in enumerate(zip(layers, lam)):
        th = L["maechtigkeit_m"]
        ax.add_patch(Rectangle((X0, d), X1 - X0, th, facecolor=SEQ(norm(lm)),
                               edgecolor=SURF, lw=0.9))
        ax.plot([X1 + 0.005, XT - 0.01], [d + th / 2, rows[i]],
                color=AXIS, lw=0.6, zorder=1)
        d += th
    for x, tt, ha in [(XT, "Schicht", "left"), (XT + 0.46, "d [m]", "right"),
                      (XT + 0.56, "λ", "right"), (XT + 0.66, "ρc", "right"),
                      (XT + 0.75, "n", "right")]:
        ax.text(x, -6.0, tt, fontsize=8.5, fontweight="bold", color=INK,
                ha=ha, va="bottom")
    for i, (L, e, lm) in enumerate(zip(layers, eff, lam)):
        ax.text(XT, rows[i], L["name"], fontsize=7.6, color=INK, va="center")
        for x, v, nd in [(XT + 0.46, L["maechtigkeit_m"], 0), (XT + 0.56, lm, 2),
                         (XT + 0.66, e["rho_c"] / 1e6, 2),
                         (XT + 0.75, L["porositaet"], 3)]:
            ax.text(x, rows[i], de(v, nd), fontsize=7.4, color=INK2,
                    ha="right", va="center")
    bt, bb = cfg["borehole"]["tiefe_kopf_m"], cfg["borehole"]["tiefe_fuss_m"]
    ax.add_patch(Rectangle((0.075, bt), 0.05, bb - bt, facecolor="none",
                           edgecolor=C2, lw=2.2, zorder=3))
    ax.annotate(f"Sonde\n{de(bt,0)}–{de(bb,0)} m", (0.10, (bt + bb) / 2),
                color=C2, fontsize=8.5, fontweight="bold", ha="center",
                va="center", zorder=4, bbox=dict(fc=SURF, ec="none", pad=1.5))
    for y, va, dy in [(0, "bottom", 7), (z_top, "top", -7)]:
        ax.annotate(f"Dirichlet T = {de(cfg['initial']['T_C'],0)} °C",
                    (X1 / 2, y), color=C4, fontsize=9, fontweight="bold",
                    ha="center", va=va, xytext=(0, dy),
                    textcoords="offset points")
    ax.annotate("Seitenränder ohne\nRandbedingung\n= Nullfluss",
                (X0 - 0.015, z_top * 0.62), color=INK2, fontsize=8,
                ha="right", va="center")
    ax.set_xlim(-0.22, XT + 0.82); ax.set_ylim(z_top + 10, -14)
    ax.set_ylabel("Tiefe unter Gelände [m]"); ax.set_xticks([])
    ax.yaxis.set_major_formatter(F)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.spines["left"].set_visible(True); ax.spines["left"].set_color(AXIS)
    ax.set_title("Untergrund, effektive Stoffwerte und Randbedingungen\n"
                 "λ in W/(m·K), ρc in MJ/(m³·K), n Porosität", loc="left")
    fig.tight_layout(); fig.savefig(FIG / "2_untergrund.png", dpi=150)
    plt.close(fig); print("-> 2_untergrund.png")

    # --- 3: Lastprofil ---
    lp = lastprofil(cfg); P = lp["P_kW"]
    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.2),
                           gridspec_kw={"width_ratios": [1.25, 1.0]})
    a = ax[0]; xx = np.arange(12)
    a.bar(xx, P, 0.62, color=[C1 if p > 0 else C4 for p in P])
    a.axhline(0, color=AXIS, lw=1.0)
    a.set_xticks(xx); a.set_xticklabels(MONATE, fontsize=7.5)
    a.set_ylabel("Feldleistung [kW]"); a.set_title("Monatsprofil des Feldes")
    a.set_ylim(P.min() * 1.35, P.max() * 1.55)
    a.annotate(f"Beladung  {de(lp['beladung_MWh'],0)} MWh/a",
               (0.4, P.max() * 1.25), color=C1, fontsize=8.5,
               fontweight="bold", ha="left")
    a.annotate(f"Entnahme  {de(lp['entnahme_MWh'],0)} MWh/a",
               (4.4, P.min() * 1.12), color=C4, fontsize=8.5,
               fontweight="bold", ha="left")
    a.grid(axis="y", alpha=0.9); a.set_axisbelow(True)
    b = ax[1]
    if "solar_MWh" in lp:
        vals = [lp["solar_MWh"], lp["bedarf_MWh"]]
        b.bar([0, 1], vals, 0.5, color=[C1, C4])
        for i, v in enumerate(vals):
            b.annotate(de(v, 0), (i, v), ha="center", va="bottom", fontsize=9,
                       fontweight="bold", color=[C1, C4][i], xytext=(0, 4),
                       textcoords="offset points")
        b.set_xticks([0, 1]); b.set_xticklabels(["Solarertrag", "Wärmebedarf"])
        b.set_ylabel("MWh/a"); b.set_ylim(0, max(vals) * 1.25)
        b.set_title(f"Jahresbilanz {de(lp['bilanz_MWh'],0)} MWh/a — Ausgleich "
                    f"bei Faktor {de(lp['faktor_bilanz_null'],4)}")
    b.grid(axis="y", alpha=0.9); b.set_axisbelow(True)
    fig.tight_layout(); fig.savefig(FIG / "3_lastprofil.png", dpi=150)
    plt.close(fig); print("-> 3_lastprofil.png")


def beispiellauf(cfg) -> None:
    """Zwei Einheitszellenläufe über 30 Jahre — wie gegeben und ausgeglichen.

    Kostet zusammen rund zwei Minuten und zeigt den Unterschied zwischen
    einem Speicher und einer Auskühlung.
    """
    import copy
    import pyvista as pv
    plt = _plt(); F = _fmt(); FIG = _figdir(cfg)
    faelle = [("wie gegeben", 0.0, C4),
              ("Bedarf auf Bilanz null", None, C1)]
    kurven = []
    for name, red, col in faelle:
        c = copy.deepcopy(cfg)
        c["field"]["einheitszelle"] = True
        c["netz"]["elemente_je_box"] = 2.0
        c["last"]["reduktionsgrad_prozent"] = (
            red if red is not None
            else lastprofil(cfg)["reduktion_bilanz_null"])
        d = HERE / "out_beispiel" / _safe(name)
        d.mkdir(parents=True, exist_ok=True)
        c["ausgabe"] = dict(cfg["ausgabe"]); c["ausgabe"]["ordner"] = str(d)
        t0 = time.time()
        netz_bauen(c, d)
        nf = netz_wandeln(c, d / f"{c['ausgabe']['prefix']}.msh", d)
        ogs_starten(prj_bauen(c, d, nf, lastkurve(c)))
        print(f"    {name}: {time.time()-t0:.0f} s")

        dom = pv.read(str(d / f"{c['ausgabe']['prefix']}_domain.vtu"))
        mid = np.asarray(dom.cell_data["MaterialIDs"])
        vol = np.abs(np.asarray(dom.compute_cell_sizes(
            length=False, area=False, volume=True).cell_data["Volume"]))
        sel = mid >= len(c["layers"])
        w = vol[sel] / vol[sel].sum()
        P = leistung_je_sonde_W(c)
        L_bh = c["borehole"]["tiefe_fuss_m"] - c["borehole"]["tiefe_kopf_m"]
        korr = fluid_korrektur(c)
        t, Tb, Tf = [], [], []
        for f in sorted(glob.glob(str(d / f"{c['ausgabe']['prefix']}_ts_*.vtu")),
                        key=lambda x: float(re.search(r"_t_([\d.]+)\.vtu$",
                                                      x).group(1))):
            m = pv.read(f)
            Tc = np.asarray(m.point_data_to_cell_data().cell_data["T"])
            td = float(re.search(r"_t_([\d.]+)\.vtu$", f).group(1)) / DAY
            box = float((Tc[sel] * w).sum()) - 273.15
            t.append(td / YEAR_D); Tb.append(box)
            Tf.append(box + P[(int(round(td / MONAT_D)) - 1) % 12] / L_bh * korr)
        kurven.append((name, np.array(t), np.array(Tb), np.array(Tf), col))

    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    ax.axhline(cfg["initial"]["T_C"], color=MUTED, lw=0.8, ls=":")
    ax.axhline(0.0, color=C4, lw=1.3, ls="--")
    xr = max(k[1].max() for k in kurven)
    for name, t, Tb, Tf, col in kurven:
        j = t > t.max() - 1.0
        ax.plot(t, Tb, color=col, lw=1.0)
        ax.plot(t, Tf, color=col, lw=0.8, alpha=0.45)
        ax.annotate(f"{name}\nSondenwand", (xr, Tb[j].min()), color=col,
                    fontsize=8.5, fontweight="bold", xytext=(6, 0),
                    textcoords="offset points", va="center")
        ax.annotate("Fluid", (xr, Tf[j].min()), color=col, fontsize=8,
                    alpha=0.75, xytext=(6, 0), textcoords="offset points",
                    va="center")
    # Links am Nullpunkt läge die Beschriftung auf der absinkenden roten
    # Kurve; im rechten Rand ist auf dieser Höhe frei.
    ax.annotate("Vereisungsgrenze\ngilt für das Fluid", (xr, 0.0),
                color=C4, fontsize=8.5, fontweight="bold", xytext=(6, -4),
                textcoords="offset points", va="top")
    ax.set_xlim(0, xr * 1.34)
    ax.set_xlabel("Betriebsjahr"); ax.set_ylabel("Temperatur [°C]")
    ax.xaxis.set_major_formatter(F); ax.yaxis.set_major_formatter(F)
    ax.set_title("Beispielrechnung an der Einheitszelle — Mittensonde eines "
                 "großen Feldes", loc="left")
    ax.grid(axis="y", alpha=0.9); ax.set_axisbelow(True)
    fig.tight_layout(); fig.savefig(FIG / "4_beispiel.png", dpi=150)
    plt.close(fig); print("-> 4_beispiel.png")


# ======================================================================
#  Dokument
# ======================================================================
def pdf_bauen(cfg) -> None:
    """LOESUNG.md nach PDF: pandoc -> HTML -> headless Chrome."""
    md = HERE / "LOESUNG.md"
    if not md.exists():
        print("LOESUNG.md nicht gefunden — übersprungen.")
        return
    pandoc = shutil.which("pandoc")
    chrome = next((p for p in (
        r"C:\Program Files\Google\Chrome\Application\chrome.exe",
        r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
        os.path.expandvars(r"%LOCALAPPDATA%\Google\Chrome\Application\chrome.exe"),
        "/usr/bin/google-chrome", "/usr/bin/chromium") if Path(p).exists()), None)
    if not pandoc or not chrome:
        print("pandoc oder Chrome fehlt — PDF übersprungen.")
        return
    html, pdf = HERE / "_LOESUNG.html", HERE / "LOESUNG.pdf"
    css = HERE / "_stil.css"
    css.write_text(
        "@page{size:A4;margin:18mm 16mm}html{font-size:10.5pt}"
        "body{font-family:'Segoe UI',system-ui,sans-serif;line-height:1.45;"
        "color:#14140f;max-width:none;margin:0}"
        "h1{font-size:17pt;border-bottom:2px solid #2a78d6;padding-bottom:.22em;"
        "margin:1.4em 0 .5em;page-break-after:avoid}"
        "h2{font-size:13pt;color:#184f95;margin:1.25em 0 .4em;"
        "page-break-after:avoid}"
        "code{font-family:Consolas,monospace;font-size:.88em;background:#f2f1ec;"
        "padding:.1em .3em;border-radius:3px}"
        "pre{background:#f6f5f0;border-left:3px solid #c3c2b7;padding:.6em .8em;"
        "page-break-inside:avoid;overflow-x:auto}pre code{background:none}"
        "table{border-collapse:collapse;width:100%;margin:.8em 0;font-size:.88em;"
        "page-break-inside:avoid}"
        "th{text-align:left;border-bottom:1.5px solid #52514e;padding:.35em .5em;"
        "background:#f6f5f0}td{border-bottom:1px solid #e1e0d9;padding:.3em .5em}"
        "figure{margin:1.1em 0;page-break-inside:avoid;text-align:center}"
        "img{max-width:100%;height:auto}"
        "figcaption{font-size:.82em;color:#52514e;text-align:left}"
        "blockquote{margin:.9em 0;padding:.5em .9em;background:#fdf6ee;"
        "border-left:3px solid #eb6834;page-break-inside:avoid}",
        encoding="utf-8")
    if subprocess.call([pandoc, str(md), "-o", str(html), "--standalone",
                        "--embed-resources", "--mathjax", "--css", str(css),
                        "--toc", "--toc-depth=2", "--metadata", "lang=de",
                        "--resource-path", str(HERE)]) != 0:
        return
    subprocess.call([chrome, "--headless=new", "--disable-gpu", "--no-sandbox",
                     f"--print-to-pdf={pdf}", "--no-pdf-header-footer",
                     "--virtual-time-budget=30000", html.as_uri()],
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    css.unlink(missing_ok=True); html.unlink(missing_ok=True)
    if pdf.exists():
        print(f"-> {pdf}  ({pdf.stat().st_size/1024:.0f} kB)")


# ======================================================================
#  Ablauf
# ======================================================================
def main() -> int:
    cfg = CONFIG
    ab = cfg["ablauf"]
    out = HERE / cfg["ausgabe"]["ordner"]

    if ab.get("bilanz", True):
        bericht(cfg)

    if ab.get("netz", False):
        out.mkdir(parents=True, exist_ok=True)
        print("\n[Netz] gmsh: Grundriss und Extrusion ...")
        netz_bauen(cfg, out)
        print("[Netz] Umwandlung nach VTU ...")
        nf = netz_wandeln(cfg, out / f"{cfg['ausgabe']['prefix']}.msh", out)
        kv = lastkurve(cfg)
        prj = prj_bauen(cfg, out, nf, kv)
        dt = MONAT_D * DAY / cfg["zeit"]["schritte_je_monat"]
        print(f"[Netz] {prj}")
        print(f"       {cfg['zeit']['jahre']} Jahre, "
              f"{int(np.ceil(kv['t_end']/dt))} Zeitschritte à {dt/DAY:.3f} d")

    if ab.get("ogs", False):
        prj = out / f"{cfg['ausgabe']['prefix']}.prj"
        if not prj.exists():
            print(f"\n{prj} fehlt — ablauf.netz einschalten.", file=sys.stderr)
        else:
            print("\n[OGS] Simulation startet ...")
            rc = ogs_starten(prj)
            if rc != 0:
                print(f"OGS endete mit {rc}", file=sys.stderr)
                return rc

    if ab.get("abbildungen", False):
        print("\n[Abbildungen]")
        abbildungen(cfg)

    if ab.get("beispiellauf", False):
        print("\n[Beispiellauf an der Einheitszelle]")
        beispiellauf(cfg)

    if ab.get("auswertung", False):
        print("\n[Auswertung]")
        auswerten(cfg, out)

    if ab.get("pdf", False):
        print("\n[PDF]")
        pdf_bauen(cfg)

    return 0


if __name__ == "__main__":
    sys.exit(main())
