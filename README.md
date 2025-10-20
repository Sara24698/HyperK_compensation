# Magnetic Field Compensation Simulation for the Hyper-Kamiokande Detector

This repository contains the Python code developed to simulate the **magnetic field compensation system** for the **Hyper-Kamiokande neutrino detector**, being constructed in Japan.<br>

The program computes the **magnetic field distribution inside the detector volume** considering the **geomagnetic field** and multiple environmental and structural effects that can distort it.  
The simulation implements the **Biot–Savart law** to evaluate how compensation coils correct for these distortions, aiming to minimize the residual magnetic field at the positions of the photomultiplier tubes (PMTs).

---

## 🧭 Overview

The simulation evaluates the local magnetic field at PMT positions under different conditions and coil configurations.  
It accounts for the following effects:

- 🌍 **Geomagnetic field** — Base magnetic field of the Earth at the Hyper-Kamiokande site.  
- 🧲 **Rebar (steel reinforcement)** — Local magnetic field perturbations induced by the steel structure.  
- 🧱 **Sagging effect** — Deformation-induced field deviations.  
- 🌀 **Gondola and calibration holes** — Field perturbations caused by non-magnetic voids or cavities.  
- 🧭 **Geomagnetic misalignment** — Possible deviation between the measured and true geomagnetic field direction.  

The simulation combines these components to estimate the **compensated field** produced by the coil system designed to neutralize all unwanted contributions.

---

## ⚙️ Structure

The code consists of several main components:

| File | Description |
|------|--------------|
| `HyperK_compensation.py` | Main execution script controlling all simulation parameters, computation, and data export. |
| `biotsavart.py` | Computes the magnetic field vector at any point using the **Biot–Savart law**. |
| `wire.py` | Defines the **geometry of the compensation coils** (rectangular, circular, or elliptical). |
| `efecto_pos2.csv` | Contains the precomputed **rebar field map** used when `--mode rebar` is selected. |
| `requirements.txt` | Lists all Python dependencies required to run the simulation. |

---

## 📘 Dependencies

All required Python packages are specified in requirements.txt.

To install them, run:

pip install -r requirements.txt

---

## ▶️ Usage

Run the simulation from the command line using Python 3:

`python main.py [--angulos <values>] [--mode <mode>] [--export_ef_data] [--histogram] [--export_results]`

Example:

`python main.py --angulos 0 90 180 --mode rebar --export_ef_data --histogram --export_results`

This command will generate:

- Compensacion.xlsx → Summary of results (mean, std. deviation, and PMT excesses).

- Datos_angulo_*.csv → Field data for each rotation.

- Histograma_angulo_*.png → Histograms of field intensity.

- Compensated_field.db → SQLite database containing detailed 3D field information.

---

## 🧩 Description of Key Functions

| Function                   | Description                                                                                  |
| -------------------------- | -------------------------------------------------------------------------------------------- |
| `posiciones()`             | Defines PMT coordinates and base field vectors.                                              |
| `rotacion_campo()`         | Rotates the magnetic field by user-defined angles.                                           |
| `Sistema_compensacion()`   | Constructs the coil system and computes the compensated field using Biot–Savart integration. |
| `export_efficiency_data()` | Saves detailed field data to CSV and SQL.                                                    |
| `export_resultados()`      | Exports summary statistics to Excel.                                                         |
| `hist()`                   | Generates histograms of magnetic field intensity.                                            |
| `resultados()`             | Main wrapper integrating all steps and exports.                                              |

---

## 🧮 Scientific Context

This simulation was developed as part of the research work on magnetic field characterization and compensation for the Hyper-Kamiokande experiment.
It aims to ensure an optimal environment for photomultiplier operation by reducing the residual field magnitude within the detector volume.

For a complete description of the methods, physics background, and experimental validation, refer to the associated publication:

📄 Rodríguez, S. et al. (2025) —
Magnetic field compensation in the Hyper-Kamiokande neutrino detector <br>
DOI: https://doi.org/10.1140/EPJP/S13360-025-06625-1

---

## 📊 Output Summary

| Output Type     | File                      | Description                                   |
| --------------- | ------------------------- | --------------------------------------------- |
| Results summary | `Compensacion.xlsx`       | Statistical summary (mean, σ, outlier PMTs).  |
| Field maps      | `Datos_angulo_*.csv`      | Field values for each PMT and rotation angle. |
| Visualizations  | `Histograma_angulo_*.png` | Field intensity histograms.                   |
| Database        | `Compensated_field.db`    | SQLite database for structured field storage. |

---

## 🧠 Notes
Magnetic fields are calculated in Tesla, though output values are usually converted to milligauss (mG).

The code can be adapted to new coil geometries by modifying wire.py.

Visual 3D representations of the system can be enabled by setting visualize=True in the Sistema_compensacion() function.

The default PMT coordinates and configurations correspond to the Hyper-Kamiokande baseline geometry.

---

## 👩‍🔬 Author
Developed by Sara Rodríguez Cabo
Institute of Space and Astronomical Science (ICTEA) — in collaboration with the University of Tokyo / Hyper-Kamiokande Collaboration.

For inquiries or collaborations, please refer to the corresponding author in the EPJ Plus publication

---

## 📄 License

This project is distributed for research and educational purposes.
Reproduction, modification, or commercial use of this code requires prior authorization from the author.


