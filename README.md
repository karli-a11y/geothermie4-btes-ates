# geothermie4-btes-ates

OpenGeoSys-Übungen zu **Borehole** und **Aquifer Thermal Energy Storage**
(BTES, ATES) — Vorlesung *Geothermie 4*.

📖 **Komplettes Tutorial:** siehe [`HANDBUCH.pdf`](HANDBUCH.pdf)
(Theorie, Übungen, Aufgaben, Plot-Interpretation, ~22 Seiten).

## Schnellstart

Voraussetzungen: **Python 3.10–3.12**, Windows / Linux / macOS.

```bash
git clone https://github.com/karli-a11y/geothermie4-btes-ates.git
cd geothermie4-btes-ates
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

`ogs.exe` muss im `PATH` auffindbar sein (wird vom `ogs`-Wheel
automatisch in `Scripts/` installiert — Verzeichnis ggf. zu `PATH`
hinzufügen). Windows-Hinweise zur Microsoft-Store-Python-Installation
im Handbuch, Abschnitt 5.

## Übung starten

```bash
cd btes/ex1_2d            # oder: ates/ex1_2d, btes/ex2_3d, ates/ex2_3d
python btes_radial_2d.py  # entsprechend: ates_radial_2d.py / btes_3d.py / ates_3d.py
python plot_results.py
```

Parameter werden im `CONFIG`-Block am Anfang des Sim-Skripts editiert
(Material, Betrieb, Zyklen, Geometrie, Mesh). Ergebnisbilder landen in
`figures/`.

## Experimentelle Varianten (Skizzen, OGS-Validierung ausstehend)

- `btes/ex2_3d/btes_3d_bhe.py` — BTES 3D mit OGS-Modul
  **`HEAT_TRANSPORT_BHE`** statt vereinfachter Volumen-Quelle. Sonden
  als 1D-Linien­elemente eingebettet, BHE-Typ (1U/2U/CXA/CXC),
  U-Rohr-Geometrie, Refrigerant und Steuerung über `CONFIG["bhe"]`
  einstellbar.
- `ates/ex1_2d/ates_radial_2d_line.py` — ATES 2D mit Filterstrecke
  als **1D-Linie** und `NodalSourceTerm` am obersten Knoten als
  Injektions­punkt (statt 2D-Volumen­filter).

Beide sind als *Skizzen* gekennzeichnet — Mesh- und PRJ-Generierung
funktionieren, der vollständige OGS-Lauf ist noch nicht
end-to-end validiert.

## Repo-Inhalt

| Pfad                          | Inhalt                                          |
|-------------------------------|-------------------------------------------------|
| `HANDBUCH.{md,pdf}`           | Tutorial (Theorie + Übungen + Aufgaben)         |
| `btes/`, `ates/`              | Übungs-Ordner mit Sim- und Plot-Skripten        |
| `solar_to_monthly.py`         | Solar­ertrag → Monatsprofil-Helfer (typisiert)  |
| `formulas/`                   | Gerenderte LaTeX-Formelgrafiken (im Handbuch)   |
| `figures_illustrations/`      | Schemazeichnungen (im Handbuch)                 |
| `requirements.txt`            | Python-Abhängigkeiten                            |
