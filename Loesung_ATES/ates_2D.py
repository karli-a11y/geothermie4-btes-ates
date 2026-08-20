#!/usr/bin/env python3
"""ATES 2D radialsymmetrisch, OHNE regionale Grundwasserstroemung.

Die Frage: Was leistet der Aquiferspeicher am Standort im guenstigsten Fall?
Ohne Hintergrundstroemung bleibt die eingespeicherte Waerme am Brunnen; das ist
die Obergrenze dessen, was der Speicher liefern kann.

Ein axialsymmetrisches Modell KANN keine gerichtete Stroemung abbilden - das
ist hier kein Mangel, sondern genau der gewuenschte Fall. Die Gegenprobe mit
Stroemung ist ates_3D.py.

    python ates_2D.py               # 30 Betriebsjahre (rund 1 h 45 min)
    python ates_2D.py --years 5     # kurzer Durchlauf zum Ausprobieren
    python ates_2D.py --no-run      # nur Netz + .prj ansehen, nicht rechnen

Ergebnisse landen neben dieser Datei in ergebnisse_2d/:
    ergebnisse_2d/ates2d.pvd       fuer ParaView (zieht alle Zeitschritte)
    ergebnisse_2d/figures/         Abbildungen und Pruefblatt
    ergebnisse_2d/*_kennzahlen.csv Rueckgewinnungsgrad, Deckungsgrad je Jahr

ZUM BEARBEITEN gibt es genau EINEN Block: das Dict FALL weiter unten.
Netz, Zeitschritt, Loeser und Ausgabe stehen darunter, sind eingestellt und
begruendet - die muss man fuer diese Uebung nicht anfassen. Der Motor liegt
im Unterordner modell/ und wird ebenfalls nicht angefasst.
"""
from __future__ import annotations

import sys
from pathlib import Path

# Der Motor liegt in modell/ neben dieser Datei - unabhaengig davon, aus
# welchem Verzeichnis heraus das Skript gestartet wird.
sys.path.insert(0, str(Path(__file__).resolve().parent / "modell"))

import ates_radial_2d as M

# ######################################################################
#
#   H I E R   S C H R A U B E N
#
#   Standort und Lastfall. Alles andere ist fertig eingestellt.
#
# ######################################################################

FALL = {

    # ------------------------------------------------------------------
    #  Lastprofil: zwoelf Monatsleistungen in Watt
    # ------------------------------------------------------------------
    #  P > 0 = einspeichern (laden),  P < 0 = foerdern (entladen).
    #
    #  Das ist die wichtigste Eingabe der ganzen Rechnung. Die Summe der
    #  Einspeisung geteilt durch die Summe der Entnahme DECKELT den
    #  Deckungsgrad, bevor eine einzige Zelle gerechnet wird:
    #  bei diesem Profil 3061/6721 = 45.5 %. Mehr als 45.5 % des Bedarfs
    #  kann der Speicher auch bei einem Rueckgewinnungsgrad von 100 %
    #  nicht liefern. Das Pruefblatt weist die Deckelung aus.
    #
    #  Die Werte unten sind Solarueberschuss minus Waermebedarf aus der
    #  Solarthermie-Uebung (halbe Krankenhauslast).
    "monatsleistung_W": [
        -595_850.0,   # Jan
        -523_490.0,   # Feb
        -336_450.0,   # Mrz
          +4_975.0,   # Apr
        +162_190.0,   # Mai
        +343_595.0,   # Jun
        +346_995.0,   # Jul
        +282_210.0,   # Aug
         +24_135.0,   # Sep
        -281_920.0,   # Okt
        -370_970.0,   # Nov
        -447_005.0,   # Dez
    ],

    # ------------------------------------------------------------------
    #  Betrieb
    # ------------------------------------------------------------------
    "T_injektion_C":   60.0,   # womit geladen wird
    "T_aquifer_C":     10.0,   # ungestoerte Untergrundtemperatur
    "betriebsjahre":     30,   # ueber --years ueberschreibbar
    #
    #  Einen Knopf fuer die Auslegungsspreizung gibt es bewusst NICHT.
    #  Das Modell rechnet den Massenstrom im Monatsprofil-Modus immer als
    #      mdot = P_max / (c_f * (T_injektion - T_aquifer))
    #  und liest operation.mass_flow_rate_kg_s dabei gar nicht. Dahinter
    #  steckt die Annahme, dass das gefoerderte Wasser bei
    #  Aquifertemperatur zurueckgespeist wird. Eine andere Spreizung ist
    #  mit diesem Skript nicht rechenbar - ein Knopf dafuer waere tot und
    #  wuerde nur die Konsolenausgabe veraendern.

    # ------------------------------------------------------------------
    #  Aquifer  (Standort aus der Abgabe)
    # ------------------------------------------------------------------
    "aquifer": {
        "maechtigkeit_m":                 38.0,
        # Durchlaessigkeitsbeiwert, der gemessene Standortwert. Die
        # Umrechnung auf die Permeabilitaet k = kf*mu/(rho*g) macht das
        # Skript weiter unten selbst - und zwar mit der Viskositaet, die
        # das Modell bei Aquifertemperatur wirklich verwendet. Wer hier
        # k statt kf eintraegt, baut sich denselben Fehler ein, der in
        # der Vorgaengerfassung steckte (siehe Kommentar bei _k()).
        "kf_m_s":                       6.0e-4,
        "porositaet":                   0.1191,
        "dichte_korn_kg_m3":            2760.0,
        "waermekapazitaet_korn_J_kgK":   793.0,
        "waermeleitfaehigkeit_korn_W_mK": 2.28,
    },

    # ------------------------------------------------------------------
    #  Deckgestein  (oben und unten gleich aufgebaut)
    # ------------------------------------------------------------------
    "deckgestein": {
        # 120 m je Seite. Gemessen an einem Vorlauf mit 80 m: die
        # Waermeleitfront (dT > 0.1 K) stand nach 30 Jahren 71.6 m weit
        # im oberen Deckgestein - nur 8.4 m vor dem auf T_amb fixierten
        # Modellrand. Nach unten 52.9 m. Die Asymmetrie ist der Auftrieb.
        "maechtigkeit_m":               120.0,
        "permeabilitaet_m2":          2.1e-16,
        # k_vertikal / k_horizontal. ACHTUNG beim Vergleich mit Fall 2:
        # dieser Knopf existiert nur im 2D-Modell. Das 3D-Modell schreibt
        # eine skalare Permeabilitaet, sein Deckgestein ist zwingend isotrop
        # und vertikal damit zehnmal durchlaessiger als hier. Wer diese Zeile
        # ins 3D-Handout kopiert, schreibt sie in ein Dict, das nie gelesen
        # wird - kein Fehler, keine Warnung, keine Wirkung.
        "anisotropie_vertikal":           0.1,
        "porositaet":                    0.05,
        "dichte_korn_kg_m3":           2700.0,
        "waermekapazitaet_korn_J_kgK":  900.0,
        "waermeleitfaehigkeit_korn_W_mK": 2.0,
    },

    # ------------------------------------------------------------------
    #  Brunnen
    # ------------------------------------------------------------------
    "brunnenradius_m": 0.1,
    # Die Filterstrecke geht ueber die volle Aquifermaechtigkeit. Anders
    # kann das 2D-Netz es nicht: es hat nur eine Zellspalte am Brunnen.

}

# ######################################################################
#
#   A B   H I E R   I S T   A L L E S   F E R T I G   E I N G E S T E L L T
#
#   Netz, Zeitschritt, Loeser, Ausgabe. Fuer diese Uebung nichts aendern.
#
# ######################################################################

C = M.CONFIG
G = 9.81
RHO_F = C["fluid"]["rho_ref_kg_m3"]

# --- Temperaturen -----------------------------------------------------
# Die Aquifertemperatur ist zugleich die Referenztemperatur des Fluids.
# Das ist kein Zufall, sondern Absicht: OGS rechnet die Viskositaet als
# mu(T) = mu_ref * (1 + slope*(T - T_ref)). Mit T_ref = T_aquifer gilt am
# ungestoerten Aquifer exakt mu = fluid.viscosity_Pa_s, und die Umrechnung
# kf -> k weiter unten ist ohne Umweg richtig.
T_AQ_K = FALL["T_aquifer_C"] + 273.15
C["initial"]["T_K"] = T_AQ_K
C["fluid"]["T_ref_K"] = T_AQ_K
C["operation"]["T_hot_K"] = FALL["T_injektion_C"] + 273.15
# operation.T_cold_K wird hier ABSICHTLICH nicht gesetzt: kein Modell liest
# den Schluessel. Bei Foerderung gibt es ueberhaupt keine Temperatur-
# Randbedingung (nur DirichletWithinTimeInterval waehrend der Beladung), die
# Foerdertemperatur ist rein dynamisch. Ihn zu setzen wuerde eine Steuerung
# vortaeuschen, die es nicht gibt.

MU = C["fluid"]["viscosity_Pa_s"]        # bei T_ref, hier also bei T_aquifer

# --- Plausibilitaet der Eingaben --------------------------------------
# OGS rechnet die Viskositaet als GERADE: mu(T) = mu_ref*(1 + slope*(T-T_ref)).
# Eine Gerade hat eine Nullstelle. Bei slope = -1.28e-2 liegt sie
# 1/0.0128 = 78.1 K ueber der Aquifertemperatur, hier also bei 88.1 GradC.
# Darueber wird die Viskositaet NEGATIV und das Modell rechnet stillschweigend
# Unsinn - es gibt dafuer keine Pruefung in OGS und keine Zeile im Pruefblatt.
# Schon ab etwa 80 GradC ist die Gerade deutlich zu flach. Wer HT-ATES rechnen
# will, braucht ein anderes Viskositaetsmodell, nicht nur ein hoeheres T_inj.
_T_NULL = FALL["T_aquifer_C"] - 1.0 / C["fluid"]["visc_slope_1_per_K"]
if FALL["T_injektion_C"] > _T_NULL - 10.0:
    raise ValueError(
        f"T_injektion_C = {FALL['T_injektion_C']:.0f} GradC ist zu hoch fuer "
        f"das lineare Viskositaetsmodell: mu(T) wird bei "
        f"{_T_NULL:.1f} GradC null und darueber negativ.\n"
        f"Zulaessig ist hier bis {_T_NULL - 10.0:.1f} GradC. Fuer hoehere "
        f"Injektionstemperaturen muss fluid.visc_slope_1_per_K angepasst "
        f"werden (oder auf None gesetzt, dann rechnet OGS mit konstanter "
        f"Viskositaet - was den Auftrieb unterschaetzt).")


def _k(kf_m_s: float) -> float:
    """Durchlaessigkeitsbeiwert kf [m/s] -> Permeabilitaet k [m2].

        k = kf * mu / (rho_f * g)

    Das mu ist NICHT frei waehlbar - es muss dasjenige sein, mit dem OGS
    bei Aquifertemperatur rechnet, sonst hat das Modell eine andere
    Durchlaessigkeit als der Standort. Genau das war in der
    Vorgaengerfassung passiert: dort stand k = 6.12e-11 m2 mit dem
    Kommentar "kf = 6e-4 m/s", umgerechnet war aber mit mu = 1.0e-3
    (Wasser bei 20 GradC) worden, waehrend das Modell bei 10 GradC mit
    mu = 1.3e-3 rechnet. Das Modell war damit 30 % zu dicht:
    K = 6.12e-11 * 9810 / 1.3e-3 = 4.62e-4 m/s statt 6.0e-4 m/s.
    """
    return kf_m_s * MU / (RHO_F * G)


# --- Aquifer ----------------------------------------------------------
_a = FALL["aquifer"]
C["layers"]["aquifer_thickness_m"] = _a["maechtigkeit_m"]
aq = C["materials"]["aquifer"]
aq["permeability_m2"] = _k(_a["kf_m_s"])
aq["permeability_ver_m2"] = _k(_a["kf_m_s"])          # isotrop
aq["porosity"] = _a["porositaet"]
aq["rho_s_kg_m3"] = _a["dichte_korn_kg_m3"]
aq["cp_s_J_kgK"] = _a["waermekapazitaet_korn_J_kgK"]
aq["lambda_s_W_mK"] = _a["waermeleitfaehigkeit_korn_W_mK"]

# --- Deckgestein ------------------------------------------------------
_d = FALL["deckgestein"]
C["layers"]["caprock_top_thickness_m"] = _d["maechtigkeit_m"]
C["layers"]["caprock_bottom_thickness_m"] = _d["maechtigkeit_m"]
for _key in ("caprock_top", "caprock_bottom"):
    cr = C["materials"][_key]
    cr["permeability_m2"] = _d["permeabilitaet_m2"]
    cr["permeability_ver_m2"] = _d["permeabilitaet_m2"] * _d["anisotropie_vertikal"]
    cr["porosity"] = _d["porositaet"]
    cr["rho_s_kg_m3"] = _d["dichte_korn_kg_m3"]
    cr["cp_s_J_kgK"] = _d["waermekapazitaet_korn_J_kgK"]
    cr["lambda_s_W_mK"] = _d["waermeleitfaehigkeit_korn_W_mK"]

# --- Brunnen ----------------------------------------------------------
C["well"]["r_well_m"] = FALL["brunnenradius_m"]
C["well"]["screen_top_offset_m"] = 0.0        # Filter ueber die volle Hoehe
C["well"]["screen_bottom_offset_m"] = 0.0     # (2D-Netz kann es nicht anders)
C["well"]["screen_permeability_m2"] = _k(_a["kf_m_s"])
# "fixed": der Massenstrom folgt dem Lastprofil. Die Vorgabe des Modells
# waere "demand" - dabei regelt der Brunnen die Foerderrate nach, bis der
# Bedarf gedeckt ist. Das ist hier ausdruecklich NICHT gewollt: eine
# nachgeregelte Rate verhindert, dass sich der Untergrund auflaedt, und
# druckt den Rueckgewinnungsgrad von 53 % auf 13 %.
C["well"]["production_control"] = "fixed"

# --- Lastprofil und Massenstrom --------------------------------------
P = FALL["monatsleistung_W"]
C["cycles"]["monthly_power_W"] = P
C["cycles"]["monthly_T_inj_K"] = None
C["cycles"]["n_cycles"] = FALL["betriebsjahre"]
C["cycles"]["ramp_days"] = 3.0
# Nur der Vollstaendigkeit halber gesetzt: im Monatsprofil-Modus liest das
# Modell diesen Schluessel NICHT, es rechnet mdot = P_max/(c_f*dT_ref) mit
# dT_ref = T_hot_K - initial.T_K. Hier wird genau diese Formel eingetragen,
# damit die Konsolenausgabe nicht etwas anderes behauptet als die Rechnung.
DT_REF = FALL["T_injektion_C"] - FALL["T_aquifer_C"]
C["operation"]["mass_flow_rate_kg_s"] = (
    max(abs(p) for p in P) / (C["fluid"]["cp_J_kgK"] * DT_REF))

# --- Physik, die nicht zum Standort gehoert ---------------------------
C["domain"]["r_max_m"] = 1000.0
C["time"]["gravity"] = True          # Auftrieb an. Bei dieser Durchlaessigkeit
                                     # ist er der Hauptverlustpfad, nicht die
                                     # Waermeleitung: Ra = 530 gegen kritisch 40.
C["dispersion"]["alpha_L_m"] = 2.0   # identisch zu Fall 2, damit die beiden
C["dispersion"]["alpha_T_m"] = 0.2   # Faelle vergleichbar bleiben

# Thermische Ausdehnung: rho(T) = rho_ref * (1 + beta*(T - T_ref)). Das ist
# die EINZIGE Zahl, die den Auftrieb steuert - und der Auftrieb ist hier der
# Hauptverlustpfad (Ra ~ 440 gegen kritisch 40), er begruendet die 120 m
# Deckgestein. Explizit gesetzt, damit er nicht stillschweigend geerbt wird.
# Es ist eine Linearisierung um T_ref: rho(60 GradC) = 980 statt real 983,2 -
# die treibende Dichtedifferenz ist rund 20 % zu gross. Qualitativ aendert das
# nichts, quantitativ ist der Auftriebsverlust systematisch etwas ueberschaetzt.
C["fluid"]["beta_1_per_K"] = -4.0e-4

# Iterationsbudget der Picard-Schleife. Die Vorgabe des Modells ist 20 - und
# genau hier lauert eine Falle: das 2D-Modell schreibt FixedTimeStepping, kann
# den Zeitschritt bei Nichtkonvergenz also NICHT verkleinern (das 3D-Modell
# kann es). Gleichzeitig rechnet dieser Fall 30 statt 2 Betriebsjahre. Der
# ausgelieferte Parametersatz braucht hoechstens 8 Iterationen, aber jede
# Verschaerfung (hoehere Injektionstemperatur, groesseres kf, groessere
# Monatsleistungen) kann das kippen - und der Lauf stirbt dann nach Stunden.
C["solver"]["nonlinear_iter"] = 50

# Speicherkoeffizient. OGS-HT liest `storage` AUSSCHLIESSLICH aus der
# Solid-Phase - der Fluid-Wert im CONFIG ist tote Konfiguration. Hier
# explizit gesetzt, damit er nicht stillschweigend aus dem Uebungsskript
# geerbt wird:  storage = S_s / (rho g), S_s = 1e-5 1/m (gespannter
# Aquifer, Lehrbuch) -> 1.0e-9 1/Pa.
# Die Ergebnisse haengen praktisch nicht daran: die hydraulische
# Diffusivitaet K/S_s liegt zwischen 4.6 und 460 m2/s, ein Drucksignal
# durchquert 150 m also in Sekunden bis Minuten - immer weit unter dem
# Zeitschritt von einem Tag. Das Druckfeld ist quasistationaer. Belegt an
# zwei fertigen Laeufen mit Faktor 100 Unterschied im Speicher-
# koeffizienten: beide treffen die stationaere Thiem-Vorhersage auf 6 %.
C["operation"]["solid_storage_1_per_Pa"] = 1.0e-9

# --- Netz -------------------------------------------------------------
# Vierteilig radial: Brunnen | Nahfeld | Fahnenband | Fernfeld. Mit EINER
# geometrischen Progression vom Brunnen bis zum Fahnenrand wird entweder
# der Brunnen zu grob oder der Fahnenrand zu fein. Im Vorlauf (feines Band
# nur bis 70 m) war die Zelle bei r = 140 m - genau dort, wo die 1-K-Front
# steht - 14.7 m breit; die Front lag in einer einzigen Zelle.
C["mesh"].update({
    "n_r_well":     6,       # ueber den Brunnenradius (0.1 m)
    "r_near_m":    20.0,     # Nahfeld, stark gestaffelt: 0.06 -> 1.0 m
    "n_r_near":    60,
    "bias_r_near":  1.05,
    "r_fine_m":   180.0,     # Fahnenband, nahezu gleichmaessig 1.0 -> 2.0 m
    "n_r_fine":   110,       # (die 1-K-Front steht nach 30 a bei ~150 m)
    "bias_r_fine":  1.0065,
    "n_r_far":     40,       # Fernfeld 180 -> 1000 m
    "bias_r_far":   1.088,
    "n_z_aquifer": 32,                 # 38 m / 32 = 1.19 m
    "caprock_fine_depth_m": 80.0,      # feine Zone deckt die Leitfront ab
    "n_z_caprock_fine":     40,        # 2.0 m
    "n_z_caprock_coarse":   12,
})

# --- Zeit und Ausgabe -------------------------------------------------
C["time"]["dt_seconds"] = 86400.0
C["time"]["output_every_n_steps"] = 10
# Absolut neben diese Datei, nicht relativ zum Arbeitsverzeichnis:
# sonst landen die Ergebnisse dort, wo man das Skript zufaellig
# aufgerufen hat.
C["output"]["out_dir"] = str(Path(__file__).resolve().parent / "ergebnisse_2d")
C["output"]["prefix"] = "ates2d"


if __name__ == "__main__":
    ein = sum(p for p in P if p > 0)
    aus = sum(-p for p in P if p < 0)
    kf = _a["kf_m_s"]
    print("FALL 1 - ATES 2D radialsymmetrisch, OHNE Grundwasserstroemung")
    print(f"  Aquifer {_a['maechtigkeit_m']:.0f} m, kf = {kf:.2e} m/s "
          f"(k = {_k(kf):.3e} m2), Deckgestein 2 x {_d['maechtigkeit_m']:.0f} m")
    print(f"  Einspeisung {ein/1e3:7.1f} kW-Monate | "
          f"Entnahme {aus/1e3:7.1f} kW-Monate")
    print(f"  -> Deckungsgrad ist bei {100*ein/aus:.1f} % gedeckelt")
    print(f"  Referenz-Massenstrom "
          f"{C['operation']['mass_flow_rate_kg_s']:.2f} kg/s")
    print()
    sys.exit(M.main())
