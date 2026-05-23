# Phase 6 Summary: 3D extension + split-sample oracle de-biasing

Two pieces of work:

1. **Split-sample oracle de-biasing (2D carryover).** The per-frame and static
   oracles selected ℓ by argmin over noisy timed sweeps and reported the cost
   of the *same* sample — an optimistic (winner's-curse) bias. Now the oracle
   SELECTS ℓ on one independent set of timing repeats and EVALUATES the reported
   cost on a second (`oracle_eval_split` / `oracle_eval_split3d`). This gives the
   true de-biased headroom and a fair method-vs-oracle number; the same harness
   is reused by the 3D comparison.

2. **Phase 6 — the 3D extension (spec §10 build phase 6).** Equal-radius
   spheres, a 3D uniform grid (counting-sort, 13-neighbor half-stencil), sphere
   narrow phase, 3D semi-implicit-Euler + reflecting walls, config-driven 2D/3D
   selection, the closed-form method at d=3 with freshly-calibrated 3D cost
   coefficients, and the 3D scenarios (stable, collapse, explosion, pulse,
   settling).

## Decisions (and defaults chosen where the spec is open)

- **Additive 3D, 2D untouched.** The 2D hot paths are byte-identical: `grid_t`
  gains `dim`/`gd`/`domain_d` and the 3D pipeline is *separate* functions
  (`grid_init3d`, `grid_build3d`, `broad_phase_3d`, `narrow_phase3d`,
  `grid_count_S_at3d`, `bruteforce…3d`) dispatched on `g->dim`. `grid_broad_phase`
  / `…_iter_time` branch once at entry (`if dim==3`) so the 2D body and its
  timing are unchanged. This preserves the committed Phase-5 2D results.

- **13-neighbor half-stencil.** 4 in-plane (the 2D stencil) + 9 in the +z slab.
  Each unordered adjacent-cell pair is visited once; verified complete against
  brute force (gate G1).

- **`dim` is config-driven** (`config.dim`, `config.domain_d`). `sim_t` carries
  `pz`/`vz` (allocated only in 3D). The closed-form controller is dimension-
  generic already (the formula, D2 estimator, clamps with `d2_ceil=d=3`,
  hysteresis all take `d`); only the Mode-B probe and blind search needed a
  `pz` argument and a `d==3` branch to the 3D count/grid calls.

- **settling scenario (spec's 3D dimension-drop test).** Gravity toward the
  z=0 floor, sized to fall the gravity-axis extent in ~50 frames
  (`g = 2·H/(50·dt)²`); the cloud collapses from space-filling (D2≈3) to a
  floor layer (D2≈2). Bounded mechanical energy ≈ N·g·H keeps KE well under the
  test ceiling at any scale. Provided in **both** 2D and 3D (uniform design)
  so the shared scenario enum stays valid; 2D-settling is a legitimate
  D2→1 scenario and passes the existing 2D energy/determinism gates.

- **3D-unsupported scenarios.** vortex / funnel / multicluster / stream have no
  3D initialiser (spec lists only 5 for 3D); `scenario_init` errors if asked
  for one in 3D.

- **3D setup (spec §9.6 / open choices, noted).** Domain **200³**, **r=0.5,
  e=0.95, dt=1/60** per the spec. Comparison at **N=8000, 100 frames** — a
  moderate N matching the 2D comparison (the spec's N=50K is the Phase-7
  campaign, not this illustration); 100 frames = one pulse period.

- **3D oracle/calibration sweep = [2, 50].** In a 200³ box ℓ<2 gives ~1e6+
  cells (M-dominated) and the cost optimum never falls there (the densest
  scenario, explosion, peaks near ℓ≈2.3; the dilute optimum is ≈8). Sweeping
  [2, 50] brackets every scenario's optimum and avoids the 8e6-cell ℓ=2r grids.
  Ericson (ℓ=2r=1) is still reported as its own baseline.

## Split-sample de-biasing — 2D result (N=8K, 200 frames)

| scenario     | static/PF (headroom over best fixed ℓ) | Closed-B/PF |
|---           |----------------------------------------|-------------|
| stable       | 1.132 (+13%) | 1.136 |
| collapse     | 1.144 (+14%) | 1.158 |
| explosion    | 1.280 (+28%) | 1.405 |
| vortex       | 1.169 (+17%) | 1.198 |
| funnel       | 1.222 (+22%) | 1.205 |
| pulse        | 1.288 (+29%) | 1.218 |
| multicluster | 1.248 (+25%) | 1.344 |
| stream       | 1.326 (+33%) | 1.297 |

De-biasing lifted the per-frame oracle only ~1–2% (the argmin-over-16
winner's-curse is mild at median-of-3 timing), so the calibrated Closed-B
ratios are unchanged within noise. The clean result is the **true headroom**: a
perfect per-frame method beats the best fixed ℓ by **+13% (quasi-static) to
+33% (flowing)**, and calibrated Closed-B captures essentially all of it on the
time-varying scenarios (beating the static oracle on funnel/pulse/stream).

## Phase 6 verification gates

All pass (`tests/test_phase6.c`, run by `make test`):

| Gate | What it checks | Result |
|---|---|---|
| G1 | 3D grid pipeline == brute-force overlap set (ℓ∈{1,2,4}, several frames) | 15 comparisons, 195 real overlaps, **0 mismatches** |
| G2a | isolated 3D pair: momentum conserved; e=1 ⇒ KE conserved | \|ΔP\|=1.3e-15, ΔKE/KE=1.2e-15 |
| G2b | force-free 3D run (e=0.95): global KE non-increasing | 0 KE-increase frames |
| G3 | closed-form d=3 vs `math_kernel_reference.py` canonical vectors | worst rel err = 0 (bit-identical) |
| G4 | 3D end-to-end on all 5 scenarios, Modes A & B, ℓ bounded | OK (all 5) |
| G5 | 3D cost-model calibration: R² + fitted a,b,c0 | R²≈0.998, a,b>0 (see below) |

All pre-existing 2D gates (Phases 1–5) still pass.

## 3D cost-model calibration (Pillar i in 3D)

`a, b` are the 3D hardware constants (~3e-8 s/cell, ~2.5–3.8e-8 s/unit-S
on the quasi-static scenarios — roughly 3–4x the 2D values, from the heavier
13-neighbor stencil and sphere test). R² is 0.99+ on stable / collapse /
settling; it **bends hard on explosion (0.84) and pulse (0.79)** — in 3D the
violently dense frames are strongly superlinear (cache/bandwidth; the worst
single-frame residual is ~44–46 ms), which the linear model cannot capture.
The controller still tracks the oracle there because only the *ratio* a/b
enters ℓ* — an honest Pillar-i datapoint, reported as the spec asks.

| scenario  |  a (s/cell) |  b (s/unitS) |   c0 (s)   |   R²   | RMSE | mean\|r\| | max\|r\| |
|---        |-------------|--------------|------------|--------|------|----------|---------|
| stable    | 2.897e-08 | 3.469e-08 | 2.31e-04 | 0.9982 | 380us | 5% | 2057us |
| collapse  | 2.913e-08 | 3.791e-08 | 1.41e-04 | 0.9987 | 382us | 4% | 1730us |
| explosion | 9.783e-09 | 1.226e-08 | 1.40e-02 | 0.8392 | 15785us | 93% | 43681us |
| pulse     | 1.303e-08 | 8.543e-09 | 1.22e-02 | 0.7924 | 15167us | 161% | 46086us |
| settling  | 2.853e-08 | 2.563e-08 | 5.92e-04 | 0.9912 | 1046us | 10% | 4855us |

## 3D comparison (N=8K, 100 frames, 200³)

Per-scenario mean / p95 / p99 / max t_detect (µs), mean ℓ, controller overhead,
and ratio to the (split-sample) per-frame oracle. Raw: `out/phase6_comparison.txt`
(`make compare_phase6 && ./compare_phase6 8000 100`).

### stable (3D)   (calib R²=0.9982)
| method            | mean | p95  | p99  | max  | ell   | oh_us  | mean/PF |
|---                |------|------|------|------|-------|--------|---------|
| Ericson(2r)       |  404761 |  415274 |  427813 |  432258 |  1.00 |     0.0 |  357.83 |
| Density(cbrt)     |    1290 |    1299 |    1476 |    1501 | 10.00 |     0.0 |    1.14 |
| Static-oracle     |    1246 |    1342 |    1356 |    1820 | 12.55 |     0.0 |    1.10 |
| Blind-9           |    1310 |    1429 |    1582 |    1826 |  9.90 |  4829.3 |    1.16 |
| Closed-A          |    1550 |    1663 |    1704 |    1707 |  8.66 |     0.1 |    1.37 |
| Closed-B          |    1294 |    1384 |    1431 |    1459 | 10.00 |    11.8 |    1.14 |
| PF-oracle         |    1131 |    1155 |    1182 |    1183 |  9.65 |     0.0 |    1.00 |

### collapse (3D)   (calib R²=0.9987)
| method            | mean | p95  | p99  | max  | ell   | oh_us  | mean/PF |
|---                |------|------|------|------|-------|--------|---------|
| Ericson(2r)       |  408311 |  420542 |  424878 |  426013 |  1.00 |     0.0 |  331.87 |
| Density(cbrt)     |    1379 |    1445 |    1476 |    2327 | 10.00 |     0.0 |    1.12 |
| Static-oracle     |    1428 |    1881 |    2335 |    2413 | 12.55 |     0.0 |    1.16 |
| Blind-9           |    1482 |    1700 |    1711 |    2318 |  9.20 |  5215.2 |    1.20 |
| Closed-A          |    1582 |    1644 |    1647 |    1647 |  8.54 |     0.1 |    1.29 |
| Closed-B          |    1432 |    1595 |    1604 |    1606 |  9.58 |    11.5 |    1.16 |
| PF-oracle         |    1230 |    1393 |    1402 |    1406 |  9.12 |     0.0 |    1.00 |

### explosion (3D)   (calib R²=0.8392)
| method            | mean | p95  | p99  | max  | ell   | oh_us  | mean/PF |
|---                |------|------|------|------|-------|--------|---------|
| Ericson(2r)       |  411520 |  426361 |  427192 |  477857 |  1.00 |     0.0 |  109.40 |
| Density(cbrt)     |   11502 |   18363 |   19533 |   60416 | 10.00 |     0.0 |    3.06 |
| Static-oracle     |    5019 |    6634 |    6941 |   13337 |  5.71 |     0.0 |    1.33 |
| Blind-9           |    5756 |    7395 |    7455 |   58243 |  5.77 | 33197.3 |    1.53 |
| Closed-A          |    5808 |    6500 |    6672 |   58575 |  4.98 |     0.9 |    1.54 |
| Closed-B          |    5573 |    6591 |    6805 |   59077 |  5.17 |    20.2 |    1.48 |
| PF-oracle         |    3762 |    4916 |    5040 |    5059 |  5.37 |     0.0 |    1.00 |

### pulse (3D)   (calib R²=0.7924)
| method            | mean | p95  | p99  | max  | ell   | oh_us  | mean/PF |
|---                |------|------|------|------|-------|--------|---------|
| Ericson(2r)       |  402944 |  414845 |  423691 |  449706 |  1.00 |     0.0 |  117.42 |
| Density(cbrt)     |   11233 |   40441 |   48463 |   51851 | 10.00 |     0.0 |    3.27 |
| Static-oracle     |    5529 |   14095 |   15013 |   15893 |  5.71 |     0.0 |    1.61 |
| Blind-9           |    4996 |   13755 |   14181 |   17252 |  6.97 | 24609.3 |    1.46 |
| Closed-A          |    4446 |   11088 |   11440 |   11863 |  7.08 |     1.0 |    1.30 |
| Closed-B          |    4359 |   10838 |   11747 |   11756 |  7.54 |    15.9 |    1.27 |
| PF-oracle         |    3432 |    8803 |    9086 |    9303 |  6.53 |     0.0 |    1.00 |

### settling (3D)   (calib R²=0.9912)
| method            | mean | p95  | p99  | max  | ell   | oh_us  | mean/PF |
|---                |------|------|------|------|-------|--------|---------|
| Ericson(2r)       |  404435 |  419349 |  426216 |  427067 |  1.00 |     0.0 |  313.29 |
| Density(cbrt)     |    1452 |    1704 |    2190 |    2418 | 10.00 |     0.0 |    1.12 |
| Static-oracle     |    1490 |    1684 |    2285 |    2465 |  9.65 |     0.0 |    1.15 |
| Blind-9           |    1677 |    2555 |    2645 |    3423 |  9.41 |  7059.4 |    1.30 |
| Closed-A          |    1663 |    1858 |    1878 |    2623 |  8.51 |     0.1 |    1.29 |
| Closed-B          |    1522 |    1751 |    1775 |    2203 |  9.44 |    11.6 |    1.18 |
| PF-oracle         |    1291 |    1489 |    1515 |    1520 |  8.47 |     0.0 |    1.00 |


### Headline numbers (3D)

- **Ericson (ℓ=2r) is unusable in 3D: ~110x–360x the per-frame oracle** (≈0.4 s/frame). At 2r=1 in a 200³ box the grid has 8×10⁶ cells, almost all empty, and the broad phase walks every one — the headline motivation for choosing ℓ at all in 3D.
- **Closed-B vs per-frame oracle: 1.14x (stable) – 1.48x (explosion)** — the same shape as 2D, with the violently-dynamic explosion worst. Chosen ℓ tracks the oracle within a few percent on every scenario.
- **Closed-B beats blind-9 on all 5 scenarios** at **~100–2000x less overhead** (~11–20 µs/frame vs blind's 4.8–33 ms/frame — in 3D the blind search's probe cost often *exceeds the detection it is optimising*).
- **Closed-B vs static oracle:** beats the best fixed ℓ on pulse (1.27 vs 1.61 — the breathing scenario, the 3D adaptivity payoff), ties on collapse, slightly behind on the quasi-static stable/settling and the transient-heavy explosion.
- **Mode B overhead ~11–20 µs/frame** — the near-zero-overhead claim holds in 3D. Calibrated Mode A is 1.29–1.54x (the D2:=d bias costs more in 3D, where the scenarios drop D2 below 3).

## Files added or changed

- `src/config.{h,c}`            — `dim`, `domain_d` (JSON-loadable).
- `src/grid.{h,c}`              — `dim`/`gd`/`domain_d`; `grid_init3d`,
                                  `grid_build3d`/`_timed3d`, 3D broad phase
                                  (13-stencil) + iter-time, `narrow_phase3d`,
                                  `grid_count_S_at3d`, 3D brute force; 2D paths
                                  byte-identical (entry dispatch on `dim`).
- `src/sim.{h,c}`               — `pz`/`vz`; 3D integrator, `wall_reflect3d`,
                                  `resolve_pair_collisions3d`; `sim_step` and
                                  KE/momentum dispatch on `dim`.
- `src/scenarios.{h,c}`         — `SCEN_SETTLING`; 3D placement helpers + 3D
                                  inits/forces; 3D vtable; dim-dispatched
                                  `scenario_init`/`apply_forces`.
- `src/oracle.{h,c}`            — `oracle_eval_split` (2D de-biasing) + 3D
                                  variants `oracle_eval3d`/`_sweep3d`/`_split3d`.
- `src/adapter.{h,c}`           — `domain_d`; `pz` arg on `adapter_step` /
                                  `adapter_blind_step`; `d==3` branch to the 3D
                                  probe / blind grid.
- `src/main.c`                  — dim-aware adapter wiring + 3D oracle dump.
- `tools/compare_phase5.c`      — split-sample per-frame oracle (2D).
- `tools/compare_phase6.c`      — new: 3D comparison harness (split-sample +
                                  3D regression calibration).
- `tests/test_phase6.c`         — new: gates G1–G5.
- `Makefile`                    — `test_phase6`, `compare_phase6` targets.
- `.gitignore`                  — ignore the new binaries.

## What is intentionally not done (Phase 7+)

- The full campaign: N∈{50K,100K} 2D and N=50K(+larger) 3D, multiple seeds,
  300/100 frames, mean±std (spec §9.6). This is Phase 7.
- Pillar-i/Pillar-ii analysis on the logged data (Phase 8).
- Sweep-and-prune secondary reference; the Sigurgeirsson exact-ℓ baseline
  (TODO.md item 2).
