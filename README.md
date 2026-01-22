# Embedded Control Benchmark: Matrix Inversion vs. O(1) Safety Invariants

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python](https://img.shields.io/badge/python-3.8+-green.svg)
![Platform](https://img.shields.io/badge/platform-ESP32%20%2F%20Embedded-orange.svg)
![Status](https://img.shields.io/badge/status-Proof%20of%20Concept-yellow.svg)

## 📉 Executive Summary

This repository contains a hardware constraint simulation demonstrating the **"Computational Gap"** in safety-critical control systems. 

It compares two approaches to implementing rigorous System Identification (SPS/LSCR) on low-power hardware (e.g., ESP32, STM32):
1.  **Traditional Online Approach:** Running complex linear algebra (`np.linalg.inv`) directly in the control loop.
2.  **MetaSpace Approach:** "Shifting Left" the complexity by pre-compiling mathematical guarantees into **O(1)** scalar checks.

The simulation proves that while advanced control theory is mathematically sound, its direct implementation on constrained hardware often leads to **Watchdog Timer (WDT) resets and catastrophic control loss.**

---

## ⚡ The "Kolumbán Paradox"

Dr. Sándor Kolumbán's PhD thesis (*System Identification in Highly Non-informative Environment*) describes brilliant algorithms for guaranteeing system safety with minimal data. However, applying these algorithms in real-time presents a hardware paradox:

* **The Math:** Requires recursive matrix operations (Computational Complexity: $O(n^3)$).
* **The Hardware:** Standard industrial microcontrollers (e.g., ESP32) have limited FLOPs and strict real-time deadlines (e.g., 20ms control loop).

If the computation time ($T_{compute}$) exceeds the loop time ($T_{loop}$), the **Watchdog Timer** resets the CPU, causing the drone/machine to crash.

---

## 🛸 Simulation Demo

![Simulation Preview](demo.gif)

<img src="crash_demo.gif" width="100%" alt="Simulation Preview">

The included script `simulation.py` visualizes two drones flying through a turbulence zone (stress test).

### Left Drone: "Traditional Implementation"
* **Logic:** Attempts to run the full System Identification math (Matrix Inversion) online during flight.
* **Behavior:** * As turbulence hits, the matrix complexity spikes.
    * CPU Load exceeds 100%.
    * **Result:** Watchdog Reset -> Motors freeze -> **CRASH.**

### Right Drone: "MetaSpace O(1) Approach"
* **Logic:** Runs pre-calculated invariant checks (Polyhedral constraints) generated offline.
* **Behavior:**
    * Checks simple inequalities: `if state < limit`.
    * CPU Load remains negligible (<1%).
    * **Result:** Deterministic safety -> **STABLE FLIGHT.**

---

## 🔬 Deep Dive: Verification

Do you want to understand the exact math behind the crash? 
[👉 Read the Technical Validation Protocol (TECHNICAL_VALIDATION.md)](TECHNICAL_VALIDATION.md)

---

## 🛠️ Installation & Usage

### Prerequisites
* Python 3.8+
* `numpy`
* `matplotlib`

### Quick Start

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/LemonScripter/embedded-control-benchmark.git](https://github.com/LemonScripter/embedded-control-benchmark.git)
    cd embedded-control-benchmark
    ```

2.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

3.  **Run the simulation:**
    ```bash
    python simulation.py
    ```

*(A window will open showing the real-time telemetry and flight path comparison.)*

---

## 📊 Technical Metrics (Simulated on ESP32 Profile)

The simulation estimates hardware load based on standard Xtensa LX6 (ESP32) performance characteristics.

| Metric | Traditional (Online Math) | MetaSpace (O(1) Check) |
| :--- | :--- | :--- |
| **Algorithmic Complexity** | $O(n^3)$ (Matrix Inv) | $O(1)$ (Scalar Logic) |
| **FLOPs per Cycle** | ~250,000 | ~20 |
| **Loop Execution Time** | > 25ms (Variable) | < 0.01ms (Deterministic) |
| **CPU Load (240MHz)** | **120% (Overload)** | **0.1%** |
| **Watchdog Risk** | CRITICAL | ZERO |

---

## 🔬 Scientific Background & Citation

This benchmark validates the engineering necessity of decoupling **mathematical verification** from **runtime execution**. It is built upon the theoretical foundations laid out in:

> **Sándor Kolumbán (2016).** *System Identification in Highly Non-informative Environment.* PhD Thesis.
> Budapest University of Technology and Economics / Vrije Universiteit Brussel.

**Note:** The MetaSpace protocol does not "simplify" the math; it moves the heavy lifting to the **Compile Time**, ensuring that the rigorous guarantees of Dr. Kolumbán's work can be deployed on affordable, low-energy chips.

---

## ⚠️ Disclaimer

**This is a Hardware Constraint Simulation.** It serves as a Proof-of-Concept (PoC) for computational cost analysis. This repository **does not** contain the proprietary source code of the MetaSpace `.bio` compiler or the polyhedral generation engine.

---


**Maintained by:** MetaSpace Technologies R&D  


