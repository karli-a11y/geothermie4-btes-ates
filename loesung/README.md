# Lösungen

Zwei Systeme, getrennt gehalten, weil sie unterschiedliche Fragen stellen und
unterschiedliche Prozesse brauchen.

| Ordner | System | Prozess | Kernfrage |
|---|---|---|---|
| [`btes/`](btes/) | Erdsonden-Wärmespeicher | `HEAT_CONDUCTION` | Trägt das Sondenfeld die geforderte Wärme über 30 Jahre, ohne einzufrieren? |
| [`ates/`](ates/) | Aquifer-Wärmespeicher | `HT` | Wieviel der eingespeicherten Wärme steht im Winter noch am Brunnen? |

Der Unterschied im Prozess ist kein Geschmack, sondern folgt aus dem Standort:
Beim BTES ist die Péclet-Zahl über den Sondenabstand rund 0,03, Advektion also
zwei Größenordnungen unter der Leitung — das Druckfeld ist dort sogar exakt
null. Beim ATES ist die Strömung genau der Punkt: die Wärmefahne wandert bei
$k_f = 6 \cdot 10^{-4}$ m/s und $i = 0{,}015$ rund 490 m im Jahr.

## BTES

| Datei | Inhalt |
|---|---|
| `btes_loesung.py` | alles in einer Datei: Bilanz, Netz, OGS-Lauf, Auswertung, Abbildungen |
| `nachhaltigkeit.py` | Reihe von Reduktionsgraden, ohne je Grad neu zu rechnen |
| `pruefung.py` | 27 Selbstprüfungen des Modells |
| `LOESUNG.md` / `.pdf` | Begleittext mit Herleitungen und Fallstricken |

Gesteuert über den `CONFIG`-Block am Dateianfang, keine Kommandozeilenschalter:

```
python btes_loesung.py
```

## ATES

| Datei | Inhalt |
|---|---|
| `ates_2D.py` | radialsymmetrisch, **ohne** regionale Strömung — die Obergrenze dessen, was der Speicher leisten kann |
| `ates_3D.py` | mit regionaler Strömung, $i = 0{,}015$ — was am Standort davon übrig bleibt |

Gesteuert über das Dict `FALL` am Dateianfang; die beiden Skripte kennen
zusätzlich Kommandozeilenschalter:

```
python ates_2D.py             # 30 Betriebsjahre
python ates_2D.py --years 5   # kurzer Durchlauf
python ates_2D.py --no-run    # nur Netz und .prj, nicht rechnen
```

## Gemeinsame Fallstricke

**Umlaute in Pfaden.** OGS öffnet auf Windows keine Projektdatei aus einem Pfad
mit Nicht-ASCII-Zeichen — der Lauf bricht mit „File does not exist" ab, obwohl
die Datei da ist. Deshalb heißt dieser Ordner `loesung` und nicht `lösung`.

**Erst klein rechnen.** Beide Systeme haben einen billigen Vorlauf: beim BTES
die Einheitszelle (`field.einheitszelle = True`, knapp eine Minute statt
Stunden), beim ATES `--years 1`. Wer eine Parameterstudie direkt auf dem vollen
Modell startet, wartet tagelang auf eine Antwort, die vorher zu haben war.

**Ergebnisordner sind nicht im Repo.** `out/`, `ergebnisse_2d/` und
`ergebnisse_3d/` erzeugt jeder Lauf selbst; sie sind bewusst ausgeschlossen,
weil sie zweistellige Gigabyte erreichen.
