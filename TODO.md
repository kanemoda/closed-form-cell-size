# Carried TODOs

Items that span phases. Each entry: when raised, what to do, by what phase.

---

## 1. Gates vs. experiment-campaign configs   (raised Phase 3)

Phase 3 verification gates run on `std_cfg` inside `tests/test_scenarios.c`:
N=400, domain 80×80, r=0.5, dt=1/60, num_frames 30–1000 depending on gate,
restitution=0.95 (G3b/G4) or 0.99 (char_pulse).

The example configs under `configs/*.json` are different — e.g. `configs/pulse.json`
uses restitution=0.99, seed=7, num_frames=200, which is *not* what G3b/G4 exercised.
Concrete divergence: G3b reported pulse max_KE = 7.37e5 (e=0.95); configs/pulse.json
peaks at 8.16e5 at frame 40 (e=0.99). The gates passed for *their* config, not
for the example config.

Spec §9.6 says the actual experiment campaign will use
N ∈ {50K, 100K}, domain 800×600, 300 frames, restitution=0.95, dt=1/60,
multiple seeds per scenario — none of which has been exercised end-to-end.

**Action (before Phase 7 / experiment harness):**
- Build a `configs/experiment/*.json` set matching §9.6 exactly.
- Re-run all Phase 1–3 gates against those configs and confirm they still pass.
- If a gate fails at experiment scale, fix or relax with justification *before*
  the campaign launches.

---

## 2. Sigurgeirsson fixed-ℓ baseline formula   (raised Phase 4)

Spec §9.2 lists "Fixed ℓ = Sigurgeirsson analytical formula" as a comparison
baseline and explicitly notes: *"TODO: extract their exact optimal-ℓ expression
when writing."* The exact paper expression is **not** in the spec.

Phase 4 currently uses a Frenkel-and-Smit-style analytical heuristic instead:
    ℓ_analytical = (V/N)^(1/d), clamped to ≥ 2r
derived by minimising a·M + b·S under uniform density with a=b. Source:
Frenkel & Smit, *Understanding Molecular Simulation* (2002), §5.2.2 —
standard cell-list balance argument.

**Action (before paper write-up):**
- Locate the Sigurgeirsson reference and either confirm the above formula is
  equivalent (it may be, under uniform-density assumptions) or replace with the
  exact published expression.
- Re-run the static-oracle vs. analytical-baseline comparison.

---

## 3. Re-run Phase 5 / 6 comparison tables after the t_narrow leak fix   (raised Phase 8)

Phase 8 Task 0 fixed a measurement leak: `sim_step` invoked
`grid_broad_phase_iter_time` (the t_M/t_S dryrun) *between* the broad-phase and
narrow-phase end-stamps, so its O(M) cost leaked into `t_narrow` and hence into
`t_detect` for every sim-measured method. The leak is now removed (the dryrun
runs after the narrow-phase stamp).

Consequence: the committed Phase 5 (`PHASE5_SUMMARY.md`) and Phase 6
(`PHASE6_SUMMARY.md`) comparison tables are slightly inflated for the
sim-measured methods (Ericson, Density, Static-oracle, Blind, Closed-A/B) by
the dryrun's O(M) cost — larger at small ell. The per-frame oracle (measured via
`oracle_eval`, which never called the dryrun) is unaffected, so the reported
gap-to-oracle ratios are slightly *over*stated.

**Action (campaign stage, NOT now):** re-run `compare_phase5` and
`compare_phase6` with the leak-fixed `sim_step` and refresh both summary tables.
Deferred deliberately so it does not delay the Phase 8 decisive test. The
qualitative conclusions are unchanged (the leak is small and hits all
sim-measured methods).
