# Lösung ATES — zwei Fälle zum selben Standort

Zwei eigenständige Rechnungen mit demselben Untergrund und demselben
Lastprofil. Sie beantworten zwei verschiedene Fragen, und **der Unterschied
zwischen ihnen ist das Ergebnis**.

| | `ates_2D.py` | `ates_3D.py` |
|---|---|---|
| Modell | 2D radialsymmetrisch (r, z) | 3D |
| Grundwasserströmung | **keine** | **i = 0,015** |
| Frage | Was leistet der Speicher im günstigsten Fall? | Was bleibt davon am Standort übrig? |
| Dauer | 30 Betriebsjahre, rund 1 h 45 min | 2 Betriebsjahre, rund 36 min |
| Zellen | 29 376 | 69 064 |
| Ergebnisse | `ergebnisse_2d/` | `ergebnisse_3d/` |

Es sind genau zwei Dateien. Der komplette Rechenkern steckt in jeder von ihnen
mit drin; zum Bearbeiten stehen die ersten rund 130 Zeilen.

```
python ates_2D.py               # kompletter Lauf
python ates_2D.py --years 5     # kurzer Durchlauf zum Ausprobieren
python ates_2D.py --no-run      # nur Netz + .prj erzeugen und ansehen

python ates_3D.py
python ates_3D.py --years 1
```

Die Skripte lassen sich aus jedem Verzeichnis starten; die Ergebnisse landen
immer neben der jeweiligen Datei. Nach jedem Lauf entstehen automatisch
Prüfblatt, Kennzahlen und Abbildungen.

---

## Was man ändert

Beide Dateien haben oben **einen** Block zum Bearbeiten:

```
######################################################################
#   H I E R   S C H R A U B E N
######################################################################
FALL = { ... }
######################################################################
#   A B   H I E R   F E R T I G   E I N G E S T E L L T
######################################################################
```

Im `FALL`-Dict steht alles, was zum Standort und zum Lastfall gehört:

| Eintrag | Bedeutung |
|---|---|
| `monatsleistung_W` | zwölf Monatsleistungen. **P > 0 = einspeichern, P < 0 = fördern** |
| `T_injektion_C`, `T_aquifer_C` | Ladetemperatur, ungestörte Untergrundtemperatur |
| `betriebsjahre` | über `--years` überschreibbar |
| `aquifer` | Mächtigkeit, `kf_m_s`, Porosität, Korndichte, -wärmekapazität, -leitfähigkeit |
| `deckgestein` | dasselbe für die Deckschichten |
| `brunnenradius_m` | bzw. `filter_kantenlaenge_m` in 3D |
| `gw_gradient`, `gw_richtung_grad` | nur `ates_3D.py` |

Alles darunter — Netz, Zeitschritt, Löser, Dispersivität, Ausgabe — ist
eingestellt und im Skript begründet. Darunter steht dann hinter einem zweiten
Balken der komplette Rechenkern (Netzgenerator, OGS-Projektdatei, Prüfbericht).
Er steht mit in der Datei, damit es bei zwei Dateien bleibt und man nichts
danebenlegen muss — angefasst wird er nicht.

Der Durchlässigkeitsbeiwert wird als **`kf_m_s`** eingetragen, nicht als
Permeabilität: die Umrechnung `k = kf·μ/(ρg)` macht das Skript selbst, und
zwar mit dem μ, das OGS bei Aquifertemperatur wirklich verwendet. Wer stattdessen
k einträgt, riskiert genau den Fehler, der in der Vorgängerfassung steckte —
dort war mit μ = 1,0·10⁻³ (Wasser bei 20 °C) umgerechnet worden, während das
Modell bei 10 °C mit 1,3·10⁻³ rechnet; der Aquifer war dadurch 30 % zu dicht.

**Drei Knöpfe gibt es bewusst nicht**, weil sie wirkungslos wären:

- eine Auslegungsspreizung. Das Modell rechnet im Monatsprofil-Modus immer
  `ṁ = P_max/(c_f·(T_inj − T_aquifer))` und liest `mass_flow_rate_kg_s` dabei
  gar nicht.
- `T_cold_K`. Kein Modell liest den Schlüssel. Bei Förderung gibt es überhaupt
  keine Temperatur-Randbedingung — die Fördertemperatur ist rein dynamisch,
  das ist ja die gesuchte Größe.
- eine Deckgesteins-Anisotropie in 3D. Das 3D-Modell schreibt nur eine skalare
  Permeabilität. `ates_2D.py` hat den Knopf, `ates_3D.py` kann ihn nicht haben.

**Eine harte Grenze:** OGS rechnet die Viskosität als Gerade,
μ(T) = μ_ref·(1 + slope·(T − T_ref)). Eine Gerade hat eine Nullstelle, hier bei
**88,1 °C** — darüber wird μ negativ und das Modell rechnet stillschweigend
Unsinn. Beide Skripte brechen deshalb mit einer Fehlermeldung ab, wenn
`T_injektion_C` über 78 °C liegt. Für HT-ATES braucht es ein anderes
Viskositätsmodell, nicht nur ein höheres T_inj.

---

## Die wichtigste Zahl steht im Lastprofil

Vor jeder Simulation gilt schon:

```
Deckungsgrad  ≤  Summe(Einspeisung) / Summe(Entnahme)
```

Beim mitgelieferten Profil sind das 3061/6721 = **45,5 %**. Mehr als 45,5 % des
Wärmebedarfs kann der Speicher auch bei einem Rückgewinnungsgrad von 100 %
nicht liefern. Wer den Deckungsgrad erhöhen will, muss an der Kollektorfläche
ansetzen, nicht am Aquifer. Das Prüfblatt und `2_deckungsgrad.png` weisen diese
Deckelung aus.

---

## Was dabei herauskommt

Wer die Dateien unverändert laufen lässt, muss diese Zahlen reproduzieren —
sie stehen so auch im Prüfblatt und in `*_kennzahlen.csv`:

| letztes Betriebsjahr | 2D, ohne GW (Jahr 30) | 3D, mit GW (Jahr 2) |
|---|---|---|
| Rückgewinnungsgrad η | **49,1 %** | **1,0 %** |
| Deckungsgrad | **22,1 %** | **0,5 %** |
| mittlere Fördertemperatur | **21,1 °C** | **10,2 °C** |
| Wärme im Deckgestein | 39 620 GJ | 12 116 GJ |
| Fahnenreichweite | 161 m | 603 m |

**Ohne** Grundwasserströmung trägt der Speicher gut ein Fünftel des
Wärmebedarfs — bei einer Obergrenze von 45,5 %, die im Lastprofil steckt. Er
schwingt sich langsam ein: η steigt von 13,5 % im ersten Jahr bis etwa Jahr 12
auf über 45 % und pendelt sich dann ein.

**Mit** Grundwasserströmung bleibt davon nichts. Der Standort hat kf = 6·10⁻⁴ m/s
und i = 0,015:

```
v_Darcy     = kf·i     = 0,78 m/d
v_Poren     = v/n      = 6,53 m/d
v_thermisch = v/(n·R)  = 1,34 m/d = 489 m/a      mit R = 4,87
```

Bei 489 m Drift pro Jahr hat sich die Sommerladung bis zum Winter vom Brunnen
gelöst. In `7_draufsicht.png` sieht man ab Tag 365 zwei getrennte Wärmekörper
stromab treiben, während am Brunnen wieder kaltes Grundwasser steht. Die
Fördertemperatur liegt bei 10,2 °C — 0,2 K über der ungestörten
Aquifertemperatur.

Der Unterschied zwischen 22,1 % und 0,5 % ist das eigentliche Ergebnis: **an
diesem Standort entscheidet nicht die Speicherauslegung, sondern die
Grundwasserströmung.** Wer hier einen ATES plant, braucht zuerst eine
belastbare Messung des Gradienten.

---

## Was nach dem Lauf entsteht

### ParaView

`ergebnisse_*/….pvd` öffnen — **nicht** die einzelnen `.vtu`. Die `.pvd` zieht
die ganze Zeitreihe als Animation mit rein. Felder: `T` [K], `p` [Pa],
`darcy_velocity` [m/s]. Daneben liegen die Teilgebiete als eigene `.vtu`
(`…_aquifer`, `…_caprock_top`, `…_hot_well_vol`, `…_lateral_inflow`) zum
Einblenden.

> **Farbskala fixieren.** ParaView skaliert per Default je Zeitschritt neu —
> damit sieht ein Lauf, in dem der halbe Aquifer kocht, genauso aus wie ein
> gesunder. `T` fest auf 283,15–333,15 K stellen.
>
> **Der 2D-Fall ist axialsymmetrisch.** Die Geometrie liegt in der
> (x = r, y = z)-Ebene. Für einen echten Schnitt den Filter *Rotational
> Extrusion* darüberlegen. Volumina sind dort Flächen: das Ringvolumen ist
> V = 2π·r·A.
>
> **Geschwindigkeit im Aquiferkern messen, nicht an der Schichtgrenze.** OGS
> mittelt die Sekundärgröße `darcy_velocity` knotenweise über die angrenzenden
> Elemente; an der Grenze Aquifer/Deckgestein steht nur die halbe
> Geschwindigkeit. Sollwert im Aquifer: 0,78 m/d.

### Abbildungen (`ergebnisse_*/figures/`)

| Datei | Die Frage, die sie beantwortet |
|---|---|
| `0_pruefblatt.png` | **Ist der Lauf gesund?** Elf Kennzahlen gegen ihr Band, Ampel und Klartext-Diagnose. Hier zuerst hinsehen. |
| `1_brunnentemperatur.png` | Wird die Fördertemperatur berechnet oder vorgeschrieben? |
| `2_deckungsgrad.png` | Trägt der Speicher die Last — und wo ist die Obergrenze? |
| `3_monatsbilanz.png` | In welchen Monaten fällt er aus, und warum? |
| `4_machbarkeitskette.png` | Ist die Anlage baubar? T_inj → ṁ → Filter → Druck |
| `5_energiebilanz.png` | Wohin geht die Wärme? Aquifer gegen Deckgestein |
| `6_feldschnitt.png` | **2D:** T-Feld (r, z) zu sechs Zeitpunkten des letzten Jahres |
| `7_draufsicht.png` | **3D:** Draufsicht auf halber Aquiferhöhe — hier sieht man die Fahne wegfließen |
| `8_laengsschnitt.png` | **3D:** Längsschnitt durch die Brunnenachse — Auftrieb und Drift |

Alle Feldbilder haben eine **feste** Farbskala T_amb…T_inj über sämtliche
Teilbilder, blaue Kontur = 1-K-Front, weiße Kontur = 50 °C.

Daneben liegen die Bilder, die das Modell selbst erzeugt. Bei einem davon
aufpassen: **`temperaturhub_ausnutzung.png` zeigt nicht η.** Dort steht die
temperaturbasierte Größe (T̄_Förder − T_amb)/(T_inj − T_amb), im Beispiel 16 % —
der energetische Rückgewinnungsgrad η = E_aus/E_ein im Prüfblatt liegt bei
49 %. Beide sind richtig, sie messen nur Verschiedenes: die eine, wie warm das
geförderte Wasser ist, die andere, wie viel Wärme zurückkommt.

> **Die Konturen im 3D-Längsschnitt nicht überinterpretieren.** Das Deckgestein
> ist am Brunnen mit rund 3 m aufgelöst, im Fernfeld mit bis zu 22 m. Wo die
> abgedriftete Fahne liegt, verschmiert die lineare Interpolation die 1-K-Front
> über ein einziges Element nach oben — sie läuft dort optisch fast an den
> oberen Rand. Knotenbasiert gemessen steht sie deutlich darunter. Für die
> **vertikalen** Verluste ist der 2D-Fall die belastbare Rechnung.

### Zahlen

- `*_kennzahlen.csv` — je Betriebsjahr: E_ein, E_aus, η, η feldbasiert,
  Deckungsgrad, Fördertemperatur, Wärme im Deckgestein, Fahnenreichweite
- `*_pruefblatt.csv` — die Ampeltabelle mit Schwellen und Diagnosen

---

## Das Prüfblatt lesen

`FEHLER` heißt: die Zahlen daneben sind wertlos, die Diagnosezeile nennt die
Stellschraube. `WARNUNG` heißt: ansehen und einordnen.

Zwei Warnungen sind normal und **kein** Fehler:

- **`T_min unter T_amb`** — Unterschwinger an der Wärmefront, eine Eigenschaft
  linearer finiter Elemente mit konsistenter Massenmatrix; sie treten auf,
  solange Δt < h²/(6a). Energetisch belanglos.
- **`T_max − T_inj`** — das Gegenstück nach oben. Im Strömungsfall rund 1,1 K,
  gut 2 % des Temperaturhubs: der auf T_inj geklemmte Filterkörper sitzt dort
  mitten im Durchstrom, und die Galerkin-Lösung überschwingt an dieser Kante.
  Energetisch 0,0004 % des Jahreseintrags. Im 2D-Fall ohne Durchstrom steht
  dort 6·10⁻¹² K, also Maschinengenauigkeit.

Im Strömungsfall kommt eine dritte Zeile dazu, die ausdrücklich **kein** Befund
ist: **`Dirichlet-Mehreintrag`** steht dort bei rund 200 %. Der
Grundwasserstrom spült ständig durch den auf T_inj geklemmten Filterkörper und
wird dabei aufgeheizt — der Brunnen wirkt als Wärmetauscher am Strom. Deshalb
ist η im Strömungsfall wenig aussagekräftig; belastbar ist dort die
**Fördertemperatur**.

Zwei Kennzahlen stehen absichtlich doppelt: **η** gegen den nominalen Eintrag
ṁ·c_p·ΔT und **η feldbasiert** gegen den tatsächlich gemessenen Wärmeeintrag.
η feldbasiert ist die konservativere Zahl.

### Eine Randbedingung, die man leicht übersieht

Im Strömungsfall tritt Grundwasser durch die Lateralfläche des Aquifers ein.
Diese Fläche trägt den regionalen Druckgradienten — und braucht zusätzlich eine
**Temperatur**-Randbedingung, sonst ist die Temperatur am Zustrom unbestimmt.
Bei advektionsdominierter Strömung ist das schlecht gestellt: in einem Lauf
ohne diese Bedingung schwankte der Zustromrand zwischen **1,8 und 23,7 °C**,
obwohl dort 10 °C anströmen, und 82 % des unterkühlten Volumens klebte an einem
Lateralrand. Das Modell weist die Zustromfläche deshalb als eigene Randgruppe
`lateral_inflow` aus und schreibt dort T = T_amb vor. Die **Abstrom**fläche
bleibt bewusst frei — sonst könnte die Wärmefahne das Modell nicht verlassen,
und genau das ist in diesem Fall die Aussage.

---

## Wenn ein Lauf abbricht

Steht im Prüfblatt `verworfene Zeitschritte: FEHLER`, ist der Lauf mitten in
der Rechnung gestorben und die Auswertung deckt nur den gerechneten Teil ab.
Es gibt zwei Ursachen, und sie brauchen entgegengesetzte Antworten:

**Budget zu knapp.** Die Zahl der Picard-Iterationen läuft gegen die Grenze,
das Residuum fällt dabei aber monoton. Dann hilft `solver.nonlinear_iter`
erhöhen. Beide Fälle stehen deshalb schon auf 50 statt 20.

**Grenzzyklus.** Das Residuum fällt nicht, es springt — im 3D-Strömungsfall
etwa 3,5·10⁻⁴ → 1,1·10⁻¹ → 2,1·10⁻³ → … Advektion, Auftrieb und der
Viskositätskontrast μ(60 °C) ≈ 0,36·μ(10 °C) koppeln sich gegenseitig auf.
**Hier hilft ein größeres Budget nicht** — mit 50 statt 20 Iterationen bricht
der Lauf an genau derselben Stelle ab, und ein kleinerer Zeitschritt ebenso
wenig. Im mitgelieferten Fall tritt das jenseits von etwa 2,7 Betriebsjahren
auf; die vorgesehenen 2 Jahre (733 d) bleiben mit Abstand darunter.

Welcher der beiden Fälle vorliegt, steht in `ergebnisse_*/driver.log`: dort die
Residuen der letzten Zeitschritte ansehen, bevor man an einer Schraube dreht.

---

## Dateien

```
ates_2D.py    Fall ohne Grundwasserströmung — FALL-Dict oben, Rechenkern darunter
ates_3D.py    Fall mit Grundwasserströmung  — dito
README.md     diese Datei
```

Mehr braucht es nicht: beide Dateien sind für sich lauffähig. Die Ergebnisse
entstehen beim Lauf in `ergebnisse_2d/` bzw. `ergebnisse_3d/` neben der
jeweiligen Datei und liegen bewusst nicht im Repository.

Voraussetzungen: OpenGeoSys 6.5.7 im Pfad, dazu `gmsh`, `pyvista`, `numpy`,
`matplotlib`. Der Bericht lässt sich mit `ATES_REPORT=0` abschalten.
