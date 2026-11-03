# Correction notice (2026-11-03)

The first version of this repository (January 2026) contained a visual demo, `simulation.py`, with a README and a `TECHNICAL_VALIDATION.md` that described it as a hardware-constraint benchmark. On review, the demo and the documents did not support the claims made. This version replaces them. The old files remain in the git history.

## What was wrong in the previous version

* **The displayed CPU load was not measured.** For the "traditional" drone it was `iterations × 4.0` (hard-coded); for the "MetaSpace" drone it was a constant 0.5 %. It did not depend on the matrix size or on any timer.
* **FLOPs, execution time and energy were formulas** (`flops / 5_000_000`, `0.8 W × time`), not measurements. The code comment assumed 60–80 MFLOPS for the ESP32 while the formula used 5 MFLOPS. The README's figures ("~250,000 FLOPs", "> 25 ms", "120 % CPU") were not produced by the code (the code computes about 21,600–86,400 FLOPs per frame at N = 15 and a load of 32–48 %, 112–128 % during the scripted stress).
* **The "watchdog reset" was a schedule**, not a reaction to overload: after the 100th call the iteration count was raised for 41 calls, which pushed the hard-coded load above 100 %, which froze the drone for 20 calls. The measured execution time (`perf_counter`) was recorded but never used.
* **The matrix inversions did nothing.** Their results were discarded; the control force did not depend on them. The same holds for the "pre-computed gain" on the MetaSpace side.
* **The two drones ran different controllers** (different gains), so their trajectories were not comparable.
* **The MetaSpace side performed no invariant-set check**, only `abs(pos) > 5`, and the on-screen text `State : OK` was a fixed string.
* **`TECHNICAL_VALIDATION.md` described things the code did not do**: a regressor matrix with least squares, 50 inversions per cycle, a measured `elapsed_time > 20 ms` watchdog trigger, a pre-compiled polyhedral set, and "mathematical equivalence" of the two sides. It also stated that the numbers were "not hard-coded" and were "measurable facts in the code".
* **The loop period was inconsistent** (20 ms in the text, 50 ms in the animation).
* **A licence badge (MIT) was shown without a licence file.** An MIT `LICENSE` file is now included.
* The README referenced `crash_demo.gif`, which was never in the repository, and `demo.gif` (7 MB) recorded the scripted animation.

## What this version contains instead

A real benchmark of the Kolumbán-set safety gate: a 500-drone fleet simulation against plain and robust EKFs, a closed-loop flight simulation, and instruction-count measurements of the same gate on an ESP32-C3-class core in an emulator. Every quoted number is produced by a script here and is labelled measured or estimated. The honest finding is mixed: the Kolumbán hull needs no noise model and no tuning and matches a robust EKF in safety, but it does not beat it, and it costs about 750 times more operations per step.

## Removed in this version

* `simulation.py` (the scripted animation), `demo.gif`, and the old validation text. They remain in the git history.
* Do not describe the old demo as a measurement anywhere (website, pitch material, audit documents).
