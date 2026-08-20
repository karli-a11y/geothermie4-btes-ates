#!/usr/bin/env python3
"""ATES 3D MIT regionaler Grundwasserstroemung, i = 0.015.

Die Frage: Was bleibt am Standort von der Speicherwirkung uebrig, wenn die
natuerliche Grundwasserstroemung mitrechnet?

Der Standort hat eine hohe Durchlaessigkeit (kf = 6e-4 m/s) UND einen steilen
Gradienten (i = 0.015). Daraus folgt:

    v_Darcy      = kf * i          = 6.0e-4 * 0.015 = 0.78 m/d
    v_Poren      = v_Darcy / n     = 0.78 / 0.1191  = 6.53 m/d
    v_thermisch  = v_Poren / R     = 6.53 / 4.87    = 1.34 m/d = 489 m/a
                   mit R = (rho c)_bulk / (n rho_f c_f) = 4.87

Die Waermefahne wandert also rund 490 m im Jahr. Was im Sommer eingespeichert
wird, steht im Winter einen halben Kilometer weiter stromab. Das ist KEIN
Modellfehler, sondern das Standortergebnis - und genau das soll dieser Lauf
zeigen. Deshalb reichen 1-2 Betriebsjahre: es gibt nichts, was sich ueber
Jahre aufladen koennte, die Antwort steht im ersten Winter.

Nur 3D kann das: ein axialsymmetrisches Modell hat konstruktiv keine
Vorzugsrichtung. Der Vergleichsfall ohne Stroemung ist ates_2D.py.

    python ates_3D.py             # 2 Betriebsjahre (rund 36 min)
    python ates_3D.py --years 1   # 1 Jahr, halbe Rechenzeit
    python ates_3D.py --no-run    # nur Netz + .prj ansehen, nicht rechnen

Ergebnisse landen neben dieser Datei in ergebnisse_3d/:
    ergebnisse_3d/ates3d_gw.pvd    fuer ParaView (zieht alle Zeitschritte)
    ergebnisse_3d/figures/         Abbildungen, u. a. Draufsicht und
                                   Laengsschnitt zu sechs Zeitpunkten

ZUM BEARBEITEN gibt es genau EINEN Block: das Dict FALL weiter unten.
Netz, Zeitschritt, Loeser und Ausgabe stehen darunter, sind eingestellt und
begruendet. Der Motor liegt im Unterordner modell/ und wird nicht angefasst.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "modell"))

import ates_3d as M

# ######################################################################
#
#   H I E R   S C H R A U B E N
#
#   Standort und Lastfall. Alles andere ist fertig eingestellt.
#
# ######################################################################

FALL = {

    # ------------------------------------------------------------------
    #  Das Merkmal dieses Falls: die regionale Grundwasserstroemung
    # ------------------------------------------------------------------
    #  Aufgepraegt als linearer Druckgradient auf der Lateralflaeche des
    #  Aquifers:  p(x) = p_hydrostatisch - rho_f * g * i * x.
    #  Gegenprobe im Ergebnis: v_x im ungestoerten Aquifer muss kf*i sein.
    #  ACHTUNG beim Nachmessen in ParaView: an der Grenzflaeche Aquifer/
    #  Deckgestein mittelt OGS die Sekundaergroesse darcy_velocity ueber
    #  beide Materialien - dort steht nur die HALBE Geschwindigkeit. Im
    #  Aquiferkern messen.
    "gw_gradient":        0.015,   # Standortwert
    "gw_richtung_grad":     0.0,   # 0 Grad = Stroemung nach +x

    # ------------------------------------------------------------------
    #  Lastprofil: zwoelf Monatsleistungen in Watt
    # ------------------------------------------------------------------
    #  Identisch zu Fall 1, damit die beiden Faelle vergleichbar sind.
    #  P > 0 = einspeichern (laden),  P < 0 = foerdern (entladen).
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
    "T_injektion_C":   60.0,
    "T_aquifer_C":     10.0,
    "betriebsjahre":      2,   # ueber --years ueberschreibbar
    #
    #  Einen Knopf fuer die Auslegungsspreizung gibt es bewusst NICHT: das
    #  Modell rechnet im Monatsprofil-Modus immer
    #      mdot = P_max / (c_f * (T_injektion - T_aquifer))
    #  und liest operation.mass_flow_rate_kg_s dabei gar nicht.

    # ------------------------------------------------------------------
    #  Aquifer  (Standort aus der Abgabe - identisch zu Fall 1)
    # ------------------------------------------------------------------
    "aquifer": {
        "maechtigkeit_m":                 38.0,
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
        # 60 m je Seite. Die reine Leitfront 2*sqrt(a*t) betraegt nach
        # 2 Jahren nur 15 m, gemessen wurde aber 31 m - der Auftrieb
        # traegt Waerme nach oben, die Leitungsformel unterschaetzt das.
        "maechtigkeit_m":                60.0,
        "permeabilitaet_m2":          2.1e-16,
        "porositaet":                    0.05,
        "dichte_korn_kg_m3":           2700.0,
        "waermekapazitaet_korn_J_kgK":  900.0,
        "waermeleitfaehigkeit_korn_W_mK": 2.0,
    },

    # ------------------------------------------------------------------
    #  Brunnen
    # ------------------------------------------------------------------
    "filter_kantenlaenge_m": 1.0,   # quadratischer Filterkoerper im 3D-Netz
    # Die Filterstrecke geht ueber die volle Aquifermaechtigkeit.

}

# ######################################################################
#
#   A B   H I E R   I S T   A L L E S   F E R T I G   E I N G E S T E L L T
#
#   Netz, Geometrie, Zeitschritt, Loeser, Ausgabe. Nichts aendern.
#
# ######################################################################

C = M.CONFIG
G = 9.81
RHO_F = C["fluid"]["rho_ref_kg_m3"]

# --- Temperaturen -----------------------------------------------------
# Aquifertemperatur = Referenztemperatur des Fluids, damit am ungestoerten
# Aquifer exakt mu = fluid.viscosity_Pa_s gilt (OGS rechnet
# mu(T) = mu_ref * (1 + slope*(T - T_ref))). Nur dann ist die Umrechnung
# kf -> k unten richtig.
T_AQ_K = FALL["T_aquifer_C"] + 273.15
C["initial"]["T_K"] = T_AQ_K
C["fluid"]["T_ref_K"] = T_AQ_K
C["operation"]["T_hot_K"] = FALL["T_injektion_C"] + 273.15
# operation.T_cold_K wird ABSICHTLICH nicht gesetzt: kein Modell liest den
# Schluessel. Bei Foerderung gibt es keine Temperatur-Randbedingung, die
# Foerdertemperatur ist rein dynamisch.

MU = C["fluid"]["viscosity_Pa_s"]        # bei T_ref, hier also bei T_aquifer

# --- Plausibilitaet der Eingaben --------------------------------------
# mu(T) = mu_ref*(1 + slope*(T - T_ref)) ist eine GERADE und hat eine
# Nullstelle: 1/0.0128 = 78.1 K ueber der Aquifertemperatur, hier also bei
# 88.1 GradC. Darueber wird mu negativ und das Modell rechnet stillschweigend
# Unsinn. Im Stroemungsfall trifft das doppelt, weil die Mobilitaet k/mu
# direkt in die Abstromgeschwindigkeit eingeht.
_T_NULL = FALL["T_aquifer_C"] - 1.0 / C["fluid"]["visc_slope_1_per_K"]
if FALL["T_injektion_C"] > _T_NULL - 10.0:
    raise ValueError(
        f"T_injektion_C = {FALL['T_injektion_C']:.0f} GradC ist zu hoch fuer "
        f"das lineare Viskositaetsmodell: mu(T) wird bei "
        f"{_T_NULL:.1f} GradC null und darueber negativ.\n"
        f"Zulaessig ist hier bis {_T_NULL - 10.0:.1f} GradC. Fuer hoehere "
        f"Injektionstemperaturen muss fluid.visc_slope_1_per_K angepasst "
        f"werden.")


def _k(kf_m_s: float) -> float:
    """Durchlaessigkeitsbeiwert kf [m/s] -> Permeabilitaet k [m2].

        k = kf * mu / (rho_f * g)

    Das mu muss dasjenige sein, mit dem OGS bei Aquifertemperatur rechnet,
    sonst hat das Modell eine andere Durchlaessigkeit als der Standort.
    Genau das war in der Vorgaengerfassung passiert: k = 6.12e-11 m2 mit
    dem Kommentar "kf = 6e-4 m/s", umgerechnet aber mit mu = 1.0e-3
    (Wasser bei 20 GradC), waehrend das Modell bei 10 GradC mit
    mu = 1.3e-3 rechnet. Der Aquifer war damit 30 % zu dicht - und die
    Fahne wanderte im Modell 377 m/a statt der 489 m/a des Standorts.
    """
    return kf_m_s * MU / (RHO_F * G)


# --- Regionale Stroemung ---------------------------------------------
C["regional_gw"]["enable"] = True
C["regional_gw"]["gradient_m_per_m"] = FALL["gw_gradient"]
C["regional_gw"]["direction_deg"] = FALL["gw_richtung_grad"]

# --- Aquifer ----------------------------------------------------------
_a = FALL["aquifer"]
C["layers"]["aquifer_thickness_m"] = _a["maechtigkeit_m"]
aq = C["materials"]["aquifer"]
aq["permeability_m2"] = _k(_a["kf_m_s"])
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
    cr["porosity"] = _d["porositaet"]
    cr["rho_s_kg_m3"] = _d["dichte_korn_kg_m3"]
    cr["cp_s_J_kgK"] = _d["waermekapazitaet_korn_J_kgK"]
    cr["lambda_s_W_mK"] = _d["waermeleitfaehigkeit_korn_W_mK"]

# --- Geometrie: Platz stromab ----------------------------------------
# Die Fahne driftet 489 m im Jahr, in zwei Betriebsjahren also 978 m. Der
# Brunnen steht deshalb weit stromauf, damit sie im Modell bleibt und man
# ihr beim Weglaufen zusehen kann: 1300 m frei stromab, 300 m stromauf.
C["domain"]["size_x_m"] = 2000.0        # Stroemungsrichtung
C["domain"]["size_y_m"] = 400.0         # quer - die Fahne wird schmal
C["wells"]["hot_well_xy"] = (-700.0, 0.0)

# --- Brunnen ----------------------------------------------------------
_s = FALL["filter_kantenlaenge_m"]
C["wells"]["screen_dx_m"] = _s
C["wells"]["screen_dy_m"] = _s
C["wells"]["screen_top_offset_m"] = 0.0
C["wells"]["screen_bottom_offset_m"] = 0.0
C["wells"]["screen_permeability_m2"] = _k(_a["kf_m_s"])
# "fixed": der Massenstrom folgt dem Lastprofil, der Brunnen regelt die
# Foerderrate NICHT nach. Die Vorgabe des Modells waere "demand".
C["wells"]["production_control"] = "fixed"

# --- Lastprofil und Massenstrom --------------------------------------
P = FALL["monatsleistung_W"]
C["cycles"]["monthly_power_W"] = P
C["cycles"]["monthly_T_inj_K"] = None
C["cycles"]["n_cycles"] = FALL["betriebsjahre"]
C["cycles"]["ramp_days"] = 3.0
# Im Monatsprofil-Modus liest das Modell diesen Schluessel NICHT; es rechnet
# mdot = P_max/(c_f*dT_ref) mit dT_ref = T_hot_K - initial.T_K. Hier wird
# genau diese Formel eingetragen, damit Konsole und Rechnung dasselbe sagen.
DT_REF = FALL["T_injektion_C"] - FALL["T_aquifer_C"]
C["operation"]["mass_flow_rate_kg_s"] = (
    max(abs(p) for p in P) / (C["fluid"]["cp_J_kgK"] * DT_REF))

# --- Physik, die nicht zum Standort gehoert ---------------------------
C["time"]["gravity"] = True
C["dispersion"]["alpha_L_m"] = 2.0
C["dispersion"]["alpha_T_m"] = 0.2

# Thermische Ausdehnung rho(T) = rho_ref*(1 + beta*(T - T_ref)) - die einzige
# Zahl, die den Auftrieb steuert. Explizit gesetzt, damit sie nicht
# stillschweigend geerbt wird; sie begruendet die 60 m Deckgestein.
C["fluid"]["beta_1_per_K"] = -4.0e-4

# Hinweis zum Vergleich mit Fall 1: dessen Deckgestein ist vertikal
# anisotrop (k_v = 0.1*k_h). Das 3D-Modell schreibt nur eine SKALARE
# Permeabilitaet - hier ist das Deckgestein zwingend isotrop und vertikal
# damit zehnmal durchlaessiger. Ein Schluessel "permeability_ver_m2" wuerde
# im 3D-CONFIG kommentarlos verschluckt.

# Speicherkoeffizient - siehe ausfuehrliche Begruendung in
# ates_2D.py. Explizit gesetzt, damit er nicht stillschweigend
# aus dem Uebungsskript geerbt wird.
C["operation"]["solid_storage_1_per_Pa"] = 1.0e-9

# --- Netz -------------------------------------------------------------
# Fein am Brunnen, wo Einspeicherung und Rueckgewinnung stattfinden. Wie
# scharf die Fahne 500 m stromab noch ist, aendert an der Aussage nichts -
# dass sie weglaeuft, sieht man auch auf grobem Fernfeldnetz. Im
# Laengsschnitt verschmiert die lineare Interpolation die 1-K-Kontur dort
# ueber eine ganze Zelle nach oben; das ist eine Eigenschaft der
# Darstellung, nicht der Rechnung.
C["mesh"] = {
    "size_in_well_m":          1.0,
    "size_near_wells_m":       4.0,
    "size_far_m":             35.0,
    "well_size_radius_m":     35.0,
    "well_size_radius_far_m": 160.0,
}

# Der Stroemungsfall ist deutlich nichtlinearer als der Speicherfall:
# Advektion, Auftrieb und Viskositaetskontrast koppeln sich
# (mu(60 GradC) = 4.7e-4 gegen 1.3e-3 kalt, k/mu lokal also 2.8-fach).
# Mit dem Standardbudget von 20 Picard-Iterationen brach der Lauf ab.
# ACHTUNG: jenseits von etwa 2.7 Betriebsjahren hilft auch dieses Budget
# nicht mehr - dort springt das Residuum in einen Grenzzyklus statt zu
# fallen, und ein groesseres Budget aendert daran nichts. Die zwei
# vorgesehenen Betriebsjahre bleiben mit Abstand darunter.
C["solver"]["nonlinear_iter"] = 50

# --- Zeit und Ausgabe -------------------------------------------------
C["time"]["dt_seconds"] = 86400.0
C["time"]["output_every_n_steps"] = 5    # dichte Ausgabe: die Drift ist schnell
# Absolut neben diese Datei, nicht relativ zum Arbeitsverzeichnis:
# sonst landen die Ergebnisse dort, wo man das Skript zufaellig
# aufgerufen hat.
C["output"]["out_dir"] = str(Path(__file__).resolve().parent / "ergebnisse_3d")
C["output"]["prefix"] = "ates3d_gw"


if __name__ == "__main__":
    kf = _a["kf_m_s"]
    n = _a["porositaet"]
    i = FALL["gw_gradient"]
    cf = C["fluid"]["cp_J_kgK"]
    R = ((n * RHO_F * cf + (1 - n) * _a["dichte_korn_kg_m3"]
          * _a["waermekapazitaet_korn_J_kgK"]) / (n * RHO_F * cf))
    vt = kf * i / n / R
    print("FALL 2 - ATES 3D MIT regionaler Grundwasserstroemung")
    print(f"  i = {i}, kf = {kf:.2e} m/s  (k = {_k(kf):.3e} m2), R = {R:.2f}")
    print(f"  v_Darcy {kf*i*86400:.2f} m/d | v_Poren {kf*i/n*86400:.2f} m/d "
          f"| v_thermisch {vt*86400:.2f} m/d = {vt*86400*365:.0f} m/a")
    drift = vt * 86400 * 365 * C["cycles"]["n_cycles"]
    platz = C["domain"]["size_x_m"] / 2 - C["wells"]["hot_well_xy"][0]
    print(f"  In {C['cycles']['n_cycles']} Jahren driftet die Front rund "
          f"{drift:.0f} m - Platz stromab: {platz:.0f} m "
          f"({100*drift/platz:.0f} % ausgenutzt)")
    print()
    sys.exit(M.main())
