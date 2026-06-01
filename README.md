# geothermie4-btes-ates

OpenGeoSys-Übungen zu **Borehole** und **Aquifer Thermal Energy Storage**
(BTES, ATES) — Vorlesung *Geothermie 4*.

📖 **Übungsskript:** siehe [`UEBUNGSSKRIPT.pdf`](UEBUNGSSKRIPT.pdf)
(Theorie, Übungen, Aufgaben, Plot-Interpretation, ~30 Seiten).

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
im Übungsskript, Abschnitt 5.

## Übung starten

```bash
cd btes/ex1_2d            # oder: ates/ex1_2d, btes/ex2_3d, ates/ex2_3d
python btes_radial_2d.py  # entsprechend: ates_radial_2d.py / btes_3d.py / ates_3d.py
python plot_results.py
```

Parameter werden im `CONFIG`-Block am Anfang des Sim-Skripts editiert
(Material, Betrieb, Zyklen, Geometrie, Mesh). Ergebnisbilder landen in
`figures/`.

## Erweiterte Varianten

- `btes/ex2_3d/btes_3d_bhe.py` — BTES 3D mit OGS-Modul
  **`HEAT_TRANSPORT_BHE`** statt vereinfachter Volumen-Quelle. Sonden
  als 1D-Linien­elemente eingebettet, BHE-Typ (1U/2U/CXA/CXC),
  U-Rohr-Geometrie, Refrigerant und Steuerung über `CONFIG["bhe"]`
  einstellbar. **Läuft end-to-end durch** (Default: 3×3-Feld, Typ 1U,
  ein Lade-/Entlade-Zyklus). Quantitative Validierung gegen die
  OGS-Benchmarks `BHE_1U`/`BHE_2U` steht noch aus — Ergebnisse
  qualitativ interpretieren.

## Repo-Inhalt

| Pfad                          | Inhalt                                          |
|-------------------------------|-------------------------------------------------|
| `UEBUNGSSKRIPT.{md,pdf}`      | Übungsskript (Theorie + Übungen + Aufgaben)     |
| `btes/`, `ates/`              | Übungs-Ordner mit Sim- und Plot-Skripten        |
| `formulas/`                   | Gerenderte LaTeX-Formelgrafiken (im Übungsskript)   |
| `figures_illustrations/`      | Schemazeichnungen (im Übungsskript)             |
| `requirements.txt`            | Python-Abhängigkeiten                            |
