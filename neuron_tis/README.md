# neuron_tis

## Overview

`neuron_tis` is a simulation and visualization project for **temporal interference stimulation (TIS)** applied to a single neuron model. The project combines biophysical neuron simulation (NEURON + LFPy), systematic electrode configuration sweeps, database-based result storage, and an interactive web dashboard for post-analysis.

The project consists of two main components:

1. **Simulation backend (`symmetry.py`)**
   Performs extracellular stimulation simulations under different electrode orientations and stimulation parameters, and stores the results in a SQLite database.

2. **Visualization frontend (`app.py`)**
   A Dash-based web application that visualizes membrane voltage traces and 3D electrode–neuron geometry from the simulation results.

---

## Features

* Biophysical neuron simulation using NEURON and LFPy
* Extracellular stimulation with multiple electrodes
* Parametric sweep over 3D electrode orientation (θ, ρ, roll)
* Sinusoidal stimulation with frequency offset (TIS-like setup)
* Automatic storage of simulation parameters, voltages, and electrode positions
* Interactive Dash dashboard

  * Voltage–time plots
  * 3D electrode and neuron visualization
  * Slider / angle-based selection of simulation cases

---

## Project Structure

```
neuron_tis/
├── symmetry.py          # Simulation and database generation
├── app.py               # Dash visualization app
├── model/
│   └── ball_and_stick.hoc   # Neuron morphology
├── DB/
│   └── SYMMETRY.db          # Simulation results (auto-generated)
└── README.md
```

---

## Requirements

### Python Packages

* numpy
* sqlite3 (standard library)
* NEURON
* LFPy
* pandas
* dash
* dash-bootstrap-components
* plotly

> ⚠️ **NEURON and LFPy** require a properly configured scientific Python environment and are typically installed via `pip`, `conda`, or system packages.

---

## Simulation Backend (`symmetry.py`)

### Purpose

`symmetry.py` performs neuron simulations under systematically rotated electrode configurations and stores:

* Simulation parameters
* Electrode coordinates
* Soma membrane voltage over time

into a SQLite database.

### Core Workflow

1. Load neuron morphology and initialize Hodgkin–Huxley dynamics
2. Generate a 3D rotation matrix from `(theta, ro, roll)`
3. Rotate a symmetric 6-electrode configuration
4. Apply sinusoidal extracellular stimulation
5. Run NEURON simulation
6. Store results in SQLite

### Database Schema

The database (`SYMMETRY.db`) contains three tables:

#### `TEST_PARAMETER`

| Column    | Description           |
| --------- | --------------------- |
| TEST_ID   | Unique simulation ID  |
| THETA     | Rotation angle (rad)  |
| RO        | Rotation angle (rad)  |
| ROLL      | Rotation angle (rad)  |
| AMPLITUDE | Stimulation amplitude |
| FREQUENCY | Base frequency        |
| DELTA     | Frequency offset      |

#### `TEST_VOLTAGE`

| Column  | Description                |
| ------- | -------------------------- |
| TEST_ID | Simulation ID              |
| TIME    | Time (ms)                  |
| VOLTAGE | Soma membrane voltage (mV) |

#### `ELECTRODE_PARAMETER`

| Column       | Description           |
| ------------ | --------------------- |
| TEST_ID      | Simulation ID         |
| ELECTRODE_ID | Electrode index       |
| X, Y, Z      | Electrode coordinates |

### Running the Simulation

```bash
python symmetry.py
```

This will perform a full sweep over:

* roll: 0–350° (step 10°)
* theta: 0–350° (step 10°)
* ro: 0–350° (step 10°)

⚠️ This results in **a large number of simulations** and may take a long time.

---

## Visualization Frontend (`app.py`)

### Purpose

`app.py` provides an interactive web interface to explore simulation results stored in the SQLite database.

### Features

* Slider-based selection of simulation cases
* Angle-based (θ) input
* Voltage vs. time plot
* 3D visualization of electrodes and neuron geometry

### Running the Dashboard

```bash
python app.py
```

Then open your browser at:

```
http://127.0.0.1:8050/
```

### Interface Description

* **Slider / θ input**: Select simulation by rotation angle
* **Voltage plot**: Soma membrane voltage over time
* **3D plot**:

  * Red / blue markers: electrodes with different frequencies
  * Central markers: neuron soma and apical dendrite

---

## Notes and Limitations

* The simulation is computationally expensive and not parallelized
* Database size can grow quickly for full angle sweeps
* Visualization currently assumes a fixed database schema

---

## Future Improvements (Suggested)

* Parallelization (MPI / multiprocessing)
* Support for multiple neuron morphologies
* Configurable stimulation waveforms
* Export figures and data
* Integration with HPC job workflows
