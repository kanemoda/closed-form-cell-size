# Phase 7 — Stage 1: experiment configs + validated headroom

Staged campaign. Stage 1 (this document) validates the headroom *measurement*
before any campaign compute is spent; the full multi-seed campaign (Stage 2) is
deliberately **not** started — the validated headroom determines the paper's
framing.

## Task 0 — experiment configs + gates at scale

`configs/experiment/*.json` encode the §9.6 physics (`r=0.5`, `e=0.95`,
`dt=1/60`) and geometry: 8 2D scenarios at `N=50000`, `800×600`, 300 frames; 5
3D scenarios at `N=50000`, `200³`, 100 frames. The full run matrix
(`N∈{50K,100K}`, multiple seeds) is obtained by overriding `N`/`seed`; see
`configs/experiment/README.md`.

`tools/validate_experiment.c` re-checks the Phase 1–3 gate *properties* at these
scales (not the small unit-test `std_cfg`). **All 13 configs pass:**

- **grid broad+narrow == O(N²) brute force** on a mid-run frame (2D and 3D) —
  exact overlap set, no missed/spurious pairs, at N=50K (the indexing/overflow
  stress the unit gates never reach);
- **KE finite** over the whole run; **force-free scenarios KE non-increasing**;
- **bit-identical determinism** across two independent runs.

One scaling caveat surfaced and is documented: **2D `multicluster` cannot pack
N=50K** into its 4 disks (~0.87 packing fraction; rejection sampling jams), so
its config is capped at N=20000. Running it at the full campaign N needs a
larger domain — flagged for the campaign harness (TODO.md item 1).

## Task 1 — a self-validated headroom measurement

"Headroom" = how much a per-frame-optimal ℓ beats the single best fixed ℓ.

### The inconsistency (why the old number was untrustworthy)

The earlier comparison measured the two oracles with **different primitives**:
the per-frame oracle via `oracle_eval` (a fresh grid per candidate, median over
K timing repeats), the static oracle via the sim's own `sim_step` path
(single-shot, persistent grid, and — a real bug — `grid_broad_phase_iter_time`,
the t_M/t_S dryrun, is invoked *between* the broad-phase end-stamp and the
narrow-phase end-stamp, so its O(M) cost leaks into `t_narrow`). The static path
was therefore systematically inflated relative to the per-frame path, which
made even a genuinely static distribution read a spurious **+13%** headroom.

### The fix

Derive **both** oracles from **one** shared `(frame × ℓ)` cost matrix, measured
by the **same** primitive — `oracle_eval_split` (de-biased: select ℓ on one
timing sample, report the cost on an independent second). Then:

- per-frame oracle = `mean_f eval[f][ argmin_e sel[f][e] ]`
- static oracle    = `mean_f eval[f][ e* ]`, `e* = argmin_e mean_f sel[f][e]`

They share the identical measurement, so they cannot diverge from a
measurement inconsistency — only from genuine frame-to-frame distribution
change. (`tools/headroom.c`.)

### Frozen-snapshot control (the self-test)

Feed the identical mid-run snapshot as every frame: every matrix row is the same
distribution, so the two oracles MUST coincide.

| frozen control | N=8K | N=50K |
|---|---|---|
| stable | +1.2% | −1.2% |
| vortex | +0.7% | +0.1% |

The frozen headroom is now **0 ± ~1.2%, straddling zero with no systematic
sign** — the wall-clock timing-noise floor — versus the old systematic **+13%**.
The measurement is consistent. (The residual is genuine timing noise at K=2;
it is centered on 0, not biased.)

### De-biased headroom: N=8K vs N=50K (byte-identical except N)

| scenario     | headroom N=8K | headroom N=50K |
|---           |---------------|----------------|
| stable       |        +0.14% |         +0.79% |
| collapse     |        +0.51% |         +0.23% |
| explosion    |        +0.85% |         -0.67% |
| vortex       |        +0.10% |         -1.39% |
| funnel       |        +4.36% |         +7.35% |
| pulse        |        +2.42% |         +2.18% |
| multicluster |        +0.32% | (N/A: packing) |
| stream       |        -0.00% |         +1.11% |

### Findings

- **The old +13–33% headroom was a measurement artifact.** With a consistent
  primitive the true per-frame-vs-best-fixed headroom is **small**: most
  scenarios sit within ±2% (timing noise) at both N. This reconciles the mature
  harness with the early ~2% N=50K diagnostic.
- **"stable" does NOT survive.** Its headroom is +0.1% (8K) / +0.8% (50K) — i.e.
  ~0. The +13% it showed before was entirely the measurement inconsistency, not
  real inelastic-clustering structure. The static-control scenario now behaves
  like a static control.
- **Real per-frame headroom is scenario-specific and modest.** Only **funnel**
  shows a consistent, clearly-above-noise benefit (**+4.4% / +7.4%** at 8K/50K)
  — filament formation genuinely rewards per-frame ℓ. **pulse** is modest
  (+2.4% / +2.2%). explosion and vortex are within timing noise at both N
  (an earlier one-off +3.5% on explosion was itself noise — the persisted re-run
  reads +0.9% / −0.7%). Net: on 6 of 8 scenarios per-frame adaptation and the
  best fixed ℓ are indistinguishable within timing noise.

### Implication for the paper's framing

The method's headline win is **not** "per-frame adaptation beats a fixed ℓ"
(that gap is small on most 2D scenarios). It is "**finds and tracks the
near-optimal ℓ automatically and for free**" — it matches the *static oracle*
(which itself sits within a few % of the per-frame oracle) and the expensive
blind search, at ~µs/frame overhead, with no manual tuning or oracle sweep, and
crushes the universal `ℓ=2r` default (≈9–20× in 2D, ≈110–360× in 3D). Per-frame
re-tuning adds a further few % only on strongly filament-forming scenarios.

## Not done (Stage 2 — deliberately deferred)

The full multi-seed campaign (`{scenario}×{N}×{seed}`, mean±std, 300/100 frames)
is **not** run. Stage 1's validated headroom is the gate; the framing above
should be confirmed before spending campaign compute.
