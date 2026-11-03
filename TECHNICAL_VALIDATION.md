# Technical validation protocol

This document states what each result in the README rests on, how it was produced, and how it can fail. It replaces the previous validation document, whose claims did not match the code (see `CHANGES.md`).

Legend: **MEASURED** = produced by running the code in this repository; **ESTIMATED** = computed from measured counts with stated assumptions; **NOT DONE** = not examined.

## 1. Test system and safety rule

Vertical axis of a drone, `v+ = a·v + b·u − g·Δt + w`, `Δt = 0.05 s` (20 Hz), `g·Δt = 0.4905`. `b` (thrust effectiveness, nominal 0.4905) and `a` (drag, 0.95–0.99) are unknown. At step 120 of a 400-step episode the payload changes: `b` becomes nominal × (1 + δ), δ uniform in ±0.3 (declared envelope valid) or ±0.45 (envelope violated for about one third of the drones), and `a` is redrawn. The state `v` is assumed to be measured exactly; only `(a, b)` are estimated.

Safety rule: `|a·v + b·u − g·Δt| ≤ 1 m/s` for every `(a, b)` the estimator considers plausible. The gate picks the command closest to the pilot's command `u_nom` that satisfies it (grid of 121 candidate commands in the Python simulation, bisection in the Rust implementation); if none exists it applies the command with the smallest worst-case value and counts an "emergency".

A **violating step** is evaluated on the **true** drone with the deterministic part only: `|a_true·v + b_true·u − g·Δt| > 1`. The wind `w` is excluded so that an estimator is not penalised for noise it cannot influence.

Pilot: PI controller plus feed-forward with the nominal model, tracking a piecewise-constant reference in ±1.4 m/s (the limit is 1.0, so the pilot sometimes asks for more than is allowed), with dither (sd 0.2). The first Nw steps are a pre-flight identification phase and are excluded from all metrics.

## 2. Estimators

* **Plain EKF**: recursive least squares with random-walk parameters (`P0 = diag(0.02², 0.10²)`, `Q = diag(1e-7, 1e-4)`), prior = nominal ± declared spread, 95 % confidence ellipsoid (χ² with 2 degrees of freedom). *Fixed*: noise variance R = 0.05². *Adaptive*: R updated from squared innovations with forgetting factor 0.97.
* **Robust EKF**: adaptive EKF with the innovation clipped to ±2.5·√S before the update (S = innovation variance). The constant 2.5 was fixed in advance and not tuned on the results.
* **Kolumbán set**: Sign-Perturbed Sums over the last Nw samples with m = 100, q = 5 (95 %). Each window sample gets an independent, data-independent sign column (sign table is fixed ROM data), which keeps the theorem's conditions for the sliding window and allows O(m) incremental updates.
* **Hull**: for a direction `c` in parameter space the `q`-th largest support among the 99 ellipsoids; a closed-form over-approximation of the set (no membership tests). Windows without excitation (smallest eigenvalue of `ΦᵀΦ/N` below 1e-7) yield no set; the fallback is the nominal model (hull only) or the declared envelope (hull + envelope).
* **Hull + envelope**: the hull intersected with the declared box (`a ∈ [0.95, 0.99]`, `b` = nominal ± 30 %); if the intersection is empty, the data-only hull is used.
* **Envelope only** and **oracle**: as in the README.

## 3. Fleet simulation (`kolumban_flotta.py`)

* Scenarios were fixed before running: calm Laplace noise (sd 0.05); heavy-tailed (Student-t, 1.5 degrees of freedom, scale 0.03, infinite variance); bursts (Laplace sd 0.05, three 30-step bursts at six times the amplitude); mixed fleet (each drone draws its noise type from {Laplace, Student-t, burst, strong wind} and a log-normal scale, sd 0.4).
* 500 drones per scenario, all estimators see identical wind, dither and reference sequences per drone (paired design). One common configuration for every drone.
* Metrics: mean violating steps; "persistently violating drone" = more than 2 % of its steps violate (threshold fixed in advance); mean of the worst 5 % of drones; intervention `|u − u_nom|`; tracking error against `clip(v_ref, ±1)`; emergency rate.
* Statistics: McNemar exact test on the paired per-drone indicator "persistently violating". 500 drones per cell; Wilson intervals for proportions in the logs.
* The oracle still violates about 0.9 % of steps in the heavy-tailed scenario because noise has already moved the speed before the next command is chosen; this is a floor, not an error.

## 4. Closed-loop flights (`kolumban_repules.py`)

Same model, 300 episodes per cell, Laplace and strong-wind noise, windows Nw = 12 and 20, valid and violated envelope. Logs: `results/kolumban_repules_v2_N*_*.txt` (include McNemar comparisons). The result most relevant to a reader: with a violated envelope all envelope-trusting variants fail (1.8–3.0 % violating steps) while the adaptive EKF does not (0.2–0.3 %).

## 5. Embedded cost (`fedelzeti/`, `fedelzeti_c3/`, `fedelzeti_c3_meres.py`)

**Implementation.** `no_std`-compatible Rust, f32, fixed-size arrays, no heap, no recursion, bounded loops. 2×2 closed forms for the 99 ellipsoids. Sliding-window update is incremental (O(m) = 99 × 5 multiply-adds per new sample), so the per-step cost does not depend on N.

**Numerical agreement (MEASURED).** `fedelzeti/log_N*.txt` and `fedelzeti_c3/log_N*.txt`: θ̂ within 3.6e-7, W within 1.2e-7, hull support (16 directions) within 6.9e-7 absolute, gate command within 0.02 on 20 cases (the Python grid step is 0.025), and the incremental window agrees with full recomputation within 1.0e-5 relative.

**Instruction counts (MEASURED).** The ESP32-C3-class build (`riscv32imc-unknown-none-elf`, no FPU, so every float operation is a soft-float routine) is loaded into the Unicorn RISC-V emulator. The harness calls `c3_init`, `c3_push`, `c3_derive`, `c3_gate` on the real window and the real 20 gate cases and counts every executed instruction by class (alu, load, store, branch, taken branch, jump, multiply, divide). Per control step (N = 20): `push` 30,730; `derive` 557,992; gate 312,624 (cheapest) / 317,557 (median) / 2,359,406 (worst, command clamped); typical step 906,279 and worst step 2,948,128 instructions. N = 12 and 40 give the same counts within 1 %.

**Time (ESTIMATED).** Cycles per instruction are not measured. Two assumption sets are reported: optimistic (1 cycle per instruction, divide 20) and pessimistic (load 2, taken branch +2, jump 2, multiply 4, divide 35). At 160 MHz this gives 5.7–7.2 ms (typical) and 18.4–23.5 ms (worst) per step. Instruction fetch from flash cache misses, interrupts and peripheral access are not modelled. The chip's actual multiplier/divider timing and pipeline behaviour were taken from general knowledge, not from a datasheet check; treat the interval as indicative.

**Other processors (ESTIMATED).** From the operation counts of the generic Rust build (`fedelzeti/log_N*.txt`: about 12,200 typical and 30,400 worst-case floating-point operations per step; divisions and square roots counted separately), with assumed cycle costs (1.5–3 cycles per flop on Cortex-M4F/M7 and ESP32, 14–60 cycles per division or square root). These are rough and were not emulated.

**PC timing (MEASURED).** WebAssembly (V8), f32, one core of an Intel i7-8550U: 37 µs typical and 129 µs worst per step at N = 20 (`fedelzeti/log_N20.txt`).

## 6. Threats to validity

* One-axis toy model; parameters, limits, noise scales and the pilot are the authors' choices. No real flight data.
* The estimators' tuning (EKF noise, prior spread, Huber constant) follows common defaults and was not optimised; a differently tuned EKF could move the ranking.
* The EKF compared here is a parameter EKF (linear in the parameters), the setting in which the Kolumbán theorem applies. Nonlinear-parameter models are outside this benchmark.
* "Hull ⊇ confidence set" is a proof sketch with numerical checks (0 of 36,000 mismatches in a separate numerical check script that is not part of this repository), not a machine-checked proof.
* The theorem's conditions (independent symmetric noise, external regressors) hold only approximately in closed loop; coverage in closed loop was measured, not proven.
* GPS jamming/spoofing, sensor faults, actuator faults, and biased or adversarial noise are NOT DONE.
* Timing on real hardware is NOT DONE.

## 7. Reproduction

See the README. Seeds are fixed in the scripts (fleet: `810000 + 1000·len(scenario) + drone index`; closed-loop: derived from Nw and the environment name), so reruns reproduce the logs on the same library versions (NumPy/SciPy). Only wall-clock figures (ms per step in Python) vary between machines.
