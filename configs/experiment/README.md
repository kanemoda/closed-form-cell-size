# Experiment configs (spec §9.6)

Template configs for the Phase-7 campaign. Each file fixes the §9.6 physics
(`r=0.5`, `restitution=0.95`, `dt=1/60`) and geometry:

- **2D** (`2d_*.json`): `domain 800×600`, `300 frames`, the 8 scenarios.
- **3D** (`3d_*.json`): `dim=3`, `domain 200³`, `100 frames`, the 5 scenarios
  (stable, collapse, explosion, pulse, settling).

## The full §9.6 run matrix

`{scenario} × {N} × {seed}`, all deterministic:

- **N.** 2D: `N ∈ {50000, 100000}`. 3D: `N ∈ {50000, 100000}` (50K plus one
  larger). The committed files are the **N=50000** templates; the campaign
  obtains the other N by overriding the `N` field (no other change).
- **seed.** Multiple seeds per scenario for mean ± std. The committed files use
  `seed=1`; the campaign sweeps `seed ∈ {1,2,3,…}` by overriding `seed`.

Keeping N and seed as the only swept fields makes every run byte-identical
except for the intended variable — the property the Phase-7 harness relies on.

## Caveat: 2D multicluster does not reach N=50K

`multicluster` places 4 disks of radius `0.4·(0.25·min(W,H)) = 60`. At N=50K
that is ~12500 particles per disk = ~0.87 areal packing fraction; random
rejection sampling jams well below that, so placement saturates and aborts.
`2d_multicluster.json` is therefore capped at **N=20000** (~0.35 packing,
places cleanly). To run multicluster at the full campaign N the scenario needs
a larger domain (or larger cluster radius) — a geometry change deferred to the
campaign harness, flagged here per TODO.md item 1 ("fix or relax with
justification *before* the campaign launches").

## Validating the gates at experiment scale

`tools/validate_experiment.c` re-checks the Phase 1–3 gate *properties* at these
scales — grid broad-phase == brute-force overlap set (sampled frames),
KE finite/bounded, momentum finite, and run-to-run determinism — without the
small-`std_cfg` the unit gates use. Run:

```
make validate_experiment
./validate_experiment configs/experiment/*.json
```
