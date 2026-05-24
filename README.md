# Closed-Form Adaptive Cell-Size Selection for Uniform-Grid Collision Detection

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![Language: C11](https://img.shields.io/badge/language-C11-blue.svg)

A uniform grid is the workhorse broad phase for particle collision detection,
and its whole performance hinges on one knob: the **cell size ℓ**. Too small and
you pay for a sea of empty cells; too large and every cell becomes a quadratic
pile-up of pair tests. The usual answers are a fixed heuristic (`ℓ = 2r`, or
`ℓ = (V/N)^{1/d}`) — fine for a uniform cloud, badly wrong for a clustered or
collapsing one.

This project derives the **optimal cell size in closed form** as an
information-theoretic functional of the *instantaneous* particle distribution,
identifies the **Rényi-2 correlation dimension `D₂`** as the quantity that
governs it, and tracks it at runtime with **near-zero overhead** — no search, no
tuning.

> Reference implementation (C, single-threaded, deterministic) for the paper
> *Entropy-Governed Adaptive Cell-Size Selection for Uniform-Grid Collision
> Detection*. The full derivation is in [`adaptive_cell_size_spec.md`](adaptive_cell_size_spec.md).

## The idea in one screen

One frame's grid-dependent cost splits into three sources — ℓ-independent
particle insertion, per-cell work `∝ M(ℓ)` (the number of cells), and pair tests
`∝ S(ℓ) = Σ nᵢ²` (the sum of squared cell occupancies):

```
T(ℓ) = c₀ + a·M(ℓ) + b·S(ℓ)
```

`S(ℓ) = Σ nᵢ² = N²·e^{−H₂(ℓ)}` is *exactly* the collision probability of the
discretised distribution (`H₂` = Rényi-2 entropy). In a scaling regime
`S(ℓ) ∝ ℓ^{D₂}`, where `D₂` is the correlation dimension: `D₂ = d` for a
space-filling cloud, `D₂ < d` when particles collapse onto a curve, shell, sheet
or filament. Minimising `T` gives a one-shot, Newton-like update:

```
              ┌                         ┐ 1/(d + D₂)
ℓ* = ℓ_cur ·  │  (d · a·M) / (D₂ · b·S)  │
              └                         ┘
```

The controller measures `M` and `S` for free during the grid build, estimates
`D₂` from a cheap 3-point log–log probe, and resizes the grid only when the
change clears a hysteresis band.

- **Mode A** — assume `D₂ := d`. Zero extra work; exact for space-filling scenes.
- **Mode B** — estimate `D₂` at runtime. Robust everywhere; the cost coefficients
  `a, b` are calibrated once by regression of the cost model (hardware constants,
  not per-frame quantities).

## Headline results

Ratios are **mean t_detect / per-frame oracle** (1.0 = matches the
non-realisable per-frame-optimal cell size), multi-seed mean (Phase 10):

| method | 2D, N=50K | 3D, N=50K | overhead |
|---|---|---|---|
| **Closed-B (ours)** | **0.99–1.07×** | **0.99–1.07×** | ~20–150 µs/frame |
| Closed-A (ours)     | 0.99–1.09×    | 1.01–1.08×    | <1 µs/frame |
| static oracle (best fixed ℓ) | 0.98–1.11× | 1.00–1.17× | — |
| blind 9-candidate search | 0.99–1.06× | 0.99–1.17× | **~6–33 ms/frame** |
| density baseline `(V/N)^{1/d}` | 1.00–2.27× | 0.99–2.04× | — |
| Ericson `ℓ = 2r` | 1.2–3.1× | **15–29×** | — |

- **Closed-B tracks the per-frame oracle to within a few percent** on every
  scenario and scale, beating the cost-matched **blind search at ~100–1000× less
  overhead** and never needing a sweep or a tuned constant.
- The naive `ℓ = 2r` default is **catastrophic in 3D** (15–29× — a 200³ grid at
  `ℓ = 2r` is mostly empty cells); the analytic density baseline degrades on
  every clustered/dynamic scene.
- **The thesis (does measuring `D₂` matter?), with error bars:** on
  D₂-controlled static configurations spanning `D₂ = 0.5 … d`, Mode A's gap to
  the oracle grows as `D₂` falls below `d`, while Mode B stays close. In 3D the
  effect is dramatic — at `D₂ = 0.5`, Mode A is **3.5×** the oracle while Mode B
  halves it to **2.0×**. The advantage scales with the dimensional deficit
  `d − D₂` and is invariant in `N` (it is a *dimension* effect, not a *scale*
  effect). See [`PHASE8_VERDICT.md`](PHASE8_VERDICT.md) and
  [`PHASE10_SUMMARY.md`](PHASE10_SUMMARY.md).

## Quickstart

Requirements: a C11 compiler (`gcc`), `make`, and Python 3 (for the math-kernel
reference and analysis). No external libraries.

```sh
git clone https://github.com/kanemoda/closed-form-cell-size.git
cd closed-form-cell-size

make            # build the simulator + all verification-gate binaries
make test       # run the full verification suite (should print "ALL OK")
make tools      # build the comparison / ablation tools

# run the simulator on an example config (2D pulse scenario, Mode-B controller):
./simulator configs/adapter_demo.json
```

A config is JSON; `dim: 3` plus `domain_d` switches to 3D. See
[`configs/`](configs/) for examples and [`configs/experiment/`](configs/experiment/)
for the paper's run matrix (spec §9.6).

## Repository layout

```
src/                 core library (single-threaded, deterministic)
  rng.c              xoshiro256++ — same seed ⇒ bit-identical stream
  grid.c             uniform grid: counting-sort binning, broad phase
                     (4-neighbour 2D / 13-neighbour 3D half-stencil), narrow phase
  sim.c              symplectic-Euler integrator, reflecting walls, collisions (2D/3D)
  scenarios.c        deterministic scenario generators (2D + 3D)
  oracle.c           non-perturbing per-ℓ cost probe (the oracle baselines)
  adapter.c          the controller: closed-form ℓ*, D₂ estimator, cost-model
                     regression, Mode A / Mode B, blind-search baseline
  baselines.c        Ericson / density fixed-ℓ baselines, static-oracle driver
  config.c json.c log.c main.c
tests/               verification gates (run by `make test`)
tools/               comparison & ablation harnesses (see "Reproducing")
configs/             example + experiment-matrix JSON configs
reference/           math_kernel_reference.py — hardware-independent kernel oracle
adaptive_cell_size_spec.md   the full technical specification
PHASE{5..10}_*.md            development-phase reports (methods + results)
```

## Verification

`make test` runs seven gate suites covering correctness and physics in 2D and
3D — the broad+narrow phase returns *exactly* the `O(N²)` brute-force
overlap set; momentum and (elastic) kinetic energy are conserved, KE is
non-increasing for restitution ≤ 1; the closed-form kernel matches
`reference/math_kernel_reference.py` to machine precision; the controller stays
bounded and converges on every scenario; and the cost-model regression recovers
known coefficients. The reference kernel suite runs standalone:

```sh
python3 reference/math_kernel_reference.py
```

## Reproducing the results

All tools are deterministic given a seed and write a results table to stderr.

```sh
# 2D / 3D method comparison (per scenario × method t_detect ratios + overhead):
./compare_phase5 <N> <frames> <seed>      # 8 2D scenarios
./compare_phase6 <N> <frames> <seed>      # 5 3D scenarios

# D₂-controlled ablation — Mode A vs Mode B vs ground-truth correlation dimension:
./phase8_d2 <dim> <N> <seed_offset>       # dim = 2 or 3

# self-validated headroom measurement (split-sample, de-biased per-frame oracle):
./headroom <N> <frames>

# re-check the verification-gate properties at campaign scale (N=50K/100K):
./validate_experiment configs/experiment/*.json
```

The development phases build the method incrementally and report results at each
step:

| doc | what it covers |
|---|---|
| [`PHASE5_SUMMARY.md`](PHASE5_SUMMARY.md) | the closed-form method + regression-calibrated coefficients (2D) |
| [`PHASE6_SUMMARY.md`](PHASE6_SUMMARY.md) | the 3D extension (13-neighbour stencil, spheres, settling) |
| [`PHASE7_HEADROOM.md`](PHASE7_HEADROOM.md) | self-validated headroom: split-sample oracle, frozen-control self-test |
| [`PHASE8_VERDICT.md`](PHASE8_VERDICT.md) | the decisive D₂ thesis test (Mode A vs Mode B) |
| [`PHASE9_SUMMARY.md`](PHASE9_SUMMARY.md) | at-scale campaign; the D₂ benefit is a dimension effect, not a scale effect |
| [`PHASE10_SUMMARY.md`](PHASE10_SUMMARY.md) | final multi-seed results tables (mean ± std) |

## Status & limitations

Research code, single-threaded by design (the cost model is kept honest by
measuring serial wall-clock; multi-core is used only to run independent configs
as separate processes). Known boundaries, documented in [`TODO.md`](TODO.md):

- A few scenario geometries (multicluster, vortex, stream) saturate their
  placement above ~20K–50K particles in the 800×600 domain and are capped there.
- 3D at N=100K on the densest scenarios is deferred (the per-frame-oracle
  candidate set blows past the buffer at the large-ℓ end).
- GPU generalisation is an explicit non-goal.

## License

[MIT](LICENSE) © 2026 kanemoda.

## Citation

If you use this work, please cite the accompanying paper (preprint forthcoming):

```bibtex
@misc{closedformcellsize,
  title  = {Entropy-Governed Adaptive Cell-Size Selection for Uniform-Grid Collision Detection},
  author = {kanemoda},
  year   = {2026},
  note   = {https://github.com/kanemoda/closed-form-cell-size}
}
```
