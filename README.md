# Embedded Control Benchmark: a Kolumbán-set safety gate on low-cost controllers

![License: MIT](https://img.shields.io/badge/license-MIT-blue)

**Status: research prototype. Simulation and emulator results only. No hardware measurements.**
Every number below comes from a script in this repository and a log in `results/` or the `fedelzeti*/log_*.txt` files.
Each result is labelled **MEASURED** (the script measured it), **ESTIMATED** (computed from measured counts with stated assumptions) or **NOT DONE**.
Read [`TECHNICAL_VALIDATION.md`](TECHNICAL_VALIDATION.md) for the method and its limits, and [`CHANGES.md`](CHANGES.md) for the correction of the previous version of this repository.

## What this repository is

A reproducible benchmark of one design question for drone control:

> When disturbances act on a drone (wind, sudden gusts, a payload change), which estimator, placed in front of the same deterministic safety gate, keeps the drone within its safety limit with the smallest change to the pilot's command, and can it run on a cheap microcontroller?

The test system is the vertical axis of a drone, `v+ = a·v + b·u − g·Δt + w` (20 Hz). The thrust effectiveness `b` and the drag `a` are unknown and change at a payload change. The safety rule is `|v+| ≤ 1 m/s` for **every** plausible `(a, b)`. The gate modifies the pilot's command `u_nom` to the nearest command that satisfies the rule.

Estimators feeding the gate:

| Name | What it is |
|---|---|
| plain EKF, fixed noise | parameter EKF (recursive least squares with random-walk parameters); the designer supplies the noise level |
| plain EKF, adaptive | same, noise level estimated from the residuals |
| robust EKF | adaptive EKF that clips any innovation larger than 2.5 standard deviations (Huber-type) |
| Kolumbán hull + envelope | Sign-Perturbed-Sums confidence set (Kolumbán, 2016; m = 100, q = 5, 95 %), over-approximated by a closed-form "hull", intersected with the declared payload envelope |
| envelope only | the declared envelope (payload ±30 %), no data |
| oracle | the true `(a, b)` is known (not realisable; a lower bound) |

## Results at a glance

### 1. Safety across a 500-drone fleet (MEASURED, simulation)

500 simulated drones per scenario, **one common configuration for all drones** (no per-drone tuning), window of N = 20 samples, declared envelope valid. Same wind, dither and pilot reference for every estimator (paired). Violating step = on the true drone `|v+| > 1 m/s` at the next step.

| Violating steps, % | calm (Laplace) | heavy-tailed gusts | noise bursts | mixed fleet |
|---|---|---|---|---|
| plain EKF, fixed noise | 0.13 | 10.7 | 9.4 | 7.2 |
| plain EKF, adaptive | 0.12 | 6.0 | 2.4 | 2.7 |
| **robust EKF** | 0.12 | 1.04 | 1.77 | 0.85 |
| **Kolumbán hull + envelope** | 0.35 | 1.18 | 1.87 | 1.01 |
| envelope only | 0.00 | 0.88 | 1.39 | 0.67 |
| oracle | 0.00 | 0.92 | 1.48 | 0.72 |

Source: `results/kolumban_flotta_log.txt` (`python kolumban_flotta.py 500 20`).

What it shows, and does not show:
* Plain EKFs violate the limit much more often when the noise is not Gaussian (fixed noise model 5–9×, adaptive 1.3–5× the hull's rate).
* The Kolumbán hull needs no noise model and no per-drone tuning, and reaches the level of the robust EKF. **It is not better than the robust EKF**: in the heavy-tailed and mixed scenarios the robust EKF has significantly fewer persistently violating drones (paired McNemar test, 21:3 and 20:7). In calm noise the plain and robust EKFs are better (0.12 % vs 0.35 %).
* The robust EKF needs about 40 floating-point operations per step (hand count), roughly 750 times fewer than the hull gate. The hull's value is a formal 95 % coverage guarantee without a noise model, not lower violation rates.
* The envelope-only gate is the safest when the envelope is true, at the price of about 45 % more intervention (0.34 vs 0.23 mean |u − u_nom|).

### 2. Closed-loop flights and a false envelope (MEASURED, simulation)

300 episodes per cell, payload change at step 120 (`results/kolumban_repules_v2_N*_*.txt`, `python kolumban_repules.py 300 20 ervenyes|serul`). If the true payload lies **outside** the declared envelope, every variant that trusts the envelope fails (hull + envelope 1.8–2.2 %, envelope only 2.3–3.0 % violating steps) while the adaptive EKF stays at 0.2–0.3 %. A hull used **without** the envelope is weak in closed loop (4–17 %): a window without excitation (constant speed and command) yields no confidence set.

### 3. Cost on embedded targets

The gate is implemented in Rust (`fedelzeti/`, f32, fixed arrays, no heap, no recursion, bounded loops; sliding window updated incrementally). The same source is built for an ESP32-C3-class core (RISC-V RV32IMC, **no FPU**, soft-float) in `fedelzeti_c3/` and run in the [Unicorn](https://www.unicorn-engine.org/) emulator on the real test vectors, with instruction counts per instruction class.

| Quantity (N = 20; N = 12 and 40 give the same costs) | Value |
|---|---|
| Output agrees with the Python f64 reference | θ within 3.4e-7, hull support within 6.9e-7, gate command within 0.02 (**MEASURED**) |
| Code size (`.text`, with soft-float routines) | 30,530 B (**MEASURED**) |
| Static RAM | about 13 KB in this harness (filter state 8,948 B including a 3,960 B sign table that could be kept in ROM as bits, about 0.5 KB) (**MEASURED**) |
| Instructions per control step, typical / worst case | 906,279 / 2,948,128 (**MEASURED**, emulator) |
| Time per step on an ESP32-C3 at 160 MHz, typical / worst case | 5.7–7.2 ms / 18.4–23.5 ms (**ESTIMATED**: optimistic 1 cycle per instruction, pessimistic load 2 / taken branch +2 / multiply 4 / divide 35 cycles; see `TECHNICAL_VALIDATION.md`) |
| Share of a 50 ms (20 Hz) loop, worst case | about 47 % (pessimistic assumption) (**ESTIMATED**) |
| PC, one core, WebAssembly, f32, per step: typical / worst | 37 µs / 129 µs at N = 20 (`fedelzeti/log_N20.txt`) (**MEASURED**, Intel i7-8550U) |

Other processor classes were **estimated** from operation counts (about 12,200 typical and 30,400 worst-case floating-point operations per step) with assumed cycle costs: ESP32 (Xtensa, with FPU) about 0.3 / 0.9 ms, Cortex-M4F about 0.3 / 0.8 ms, Cortex-M7 about 0.06 / 0.16 ms. These could **not** be emulated here and are not measurements.

Consequences:
* 20 Hz control is feasible on an FPU-less ESP32-C3-class core; 50 Hz or faster is not, in the worst case.
* For a fleet, every drone can run its own gate independently (cost scales linearly, nothing is shared). A central computer is also possible (one PC core about 37 µs per drone-step) but that was **not measured for 500 drones together**.

## What we do not claim

* That the Kolumbán set makes the control safer than a well-designed robust EKF (it does not, in these tests).
* That it is cheaper than an EKF (it is about 750 times more expensive per step); the cost result is that it is **affordable** at 20 Hz on a low-cost controller.
* That the 95 % level is 100 %, or that anything holds outside the theorem's conditions (linear in the parameters, independent symmetric noise).
* Anything about GPS jamming, GPS spoofing, sensor faults, rotor failure, or real flight data. The state (vertical speed) is assumed to be measured exactly. Only wind-type disturbances and payload/drag changes are modelled.
* Anything about hardware timing. No microcontroller board was used.

## Reproduce

```
pip install -r requirements.txt
# fleet simulation (about 15-20 minutes per scenario on 7 cores; quick check: 8 drones)
python kolumban_flotta.py 8 20 laplace 2
python kolumban_flotta.py 500 20 laplace,vastag,roham,vegyes 7
# closed-loop flights
python kolumban_repules.py 300 20 ervenyes
# embedded cost
rustup target add riscv32imc-unknown-none-elf
cd fedelzeti_c3 && cargo build --release && cd ..
python fedelzeti_c3_meres.py 20
```

The PC timing build is `cd fedelzeti && cargo build --release --target wasm32-wasip1` and `node run_wasi.mjs <wasm> tv20.txt 20000` (Node.js WASI). Test vectors (`fedelzeti/tv*.txt`) are produced by `fedelzeti_tesztvektor.py` (environment variable `NW`).
Code comments and log headings are in Hungarian (the scenario names in the logs: `laplace` calm, `vastag` heavy-tailed, `roham` burst, `vegyes` mixed, `ervenyes` valid envelope, `serul` violated envelope).

## Layout

| Path | Content |
|---|---|
| `kolumban_flotta.py` | 500-drone fleet simulation (scenarios fixed before running) |
| `kolumban_repules.py` | closed-loop flight simulation and all estimators |
| `kolumban_dron.py`, `kolumban_dron_elozetes.py`, `kolumban_arx.py`, `kolumban_runtime_bench.py`, `kolumban_vezerles_futasido.py`, `kolumban_decision.py` | drone model, SPS set, hull, envelope; supporting modules |
| `fedelzeti/` | Rust (f32, std) gate: validation against Python, PC timing, operation counting |
| `fedelzeti_c3/` | the same algorithm as `no_std` Rust for RISC-V RV32IMC |
| `fedelzeti_c3_meres.py` | emulator harness (Unicorn): loads the ELF, runs the test vectors, counts instructions |
| `results/` | logs of the runs quoted above |

## Scope and disclaimer

* The safety gate here is **hand-written Rust** implementing the rule `|v+| ≤ 1 m/s for all plausible (a, b)`. It is **not** generated by the MetaSpace `.bio` compiler and **has no formal (Z3) proof** in this repository. The `.bio` compiler and polytope generation engine are not part of this repository.
* The "hull ⊇ confidence set" property is a proof sketch with numerical checks, not a machine-checked proof.
* One-axis toy model, parameters chosen by the authors, simulation only.

## Scientific background

Sándor Kolumbán (2016), *System Identification in Highly Non-informative Environment*, PhD thesis, Budapest University of Technology and Economics / Vrije Universiteit Brussel (Sign-Perturbed Sums, theorem 2.17). The thesis does not claim the method is unusable at run time; it calls the exploration of the full region computationally taxing. This repository measures what a closed-form over-approximation of that region costs on small controllers.

## Licence

MIT, see [`LICENSE`](LICENSE). The previous README showed an MIT badge without a licence file; the file is now present.

Maintained by MetaSpace Technologies R&D · metaspace.bio
