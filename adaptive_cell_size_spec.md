# Technical Specification
## Entropy-Governed Adaptive Cell-Size Selection for Uniform-Grid Collision Detection

**Author:** Efe Deniz Bağlar, Department of Computer Engineering, Ege University
**Document status:** v1 — implementation & writing spec. Consolidates the theory, algorithm,
and experimental design agreed in planning. Written in English because it feeds the
English arXiv paper and the implementation.

This document has three jobs:
1. The complete build spec for a from-scratch, reproducible implementation (via Claude Code).
2. The backbone of the paper's *Method* and *Experiments* sections.
3. The single reference that prevents re-deriving or re-deciding anything.

---

## 0. Goal and Framing

**Paper type.** Theory-driven empirical systems paper. The theory is the spine; the
experiments validate it. Not a heuristic-with-benchmarks paper.

**Target.** arXiv preprint; later submission possible to a graphics venue
(MIG / SCA / EG Short Papers / JCGT).

**Core claim.** The cost-optimal uniform-grid cell size for broad-phase collision
detection of equal-radius particles is determined, in closed form, by the **correlation
(Rényi-2) dimension** of the *instantaneous* particle distribution. This yields a
runtime adaptation rule with **near-zero overhead**.

**Objective being optimized.** Measured wall-clock frame time (equivalently, FPS).
The theory provides the *mechanism*; the ground truth is always measured frame time.
Two distinct notions of "optimal ℓ" are kept strictly separate throughout:
- `ℓ_emp` — the cell size that yields the lowest *measured* frame time. Ground truth.
- `ℓ_theory` — the minimizer of the entropy-based cost model. A prediction.
The paper's empirical claims are always stated relative to `ℓ_emp`. The gap
`ℓ_theory − ℓ_emp` is a *reported result*, never hidden.

**Relation to the previous draft.** The earlier paper used a heuristic search over
7–9 candidate cell sizes, had no closed-form theory, a single oracle baseline, and
no error bars. The source code is lost; all prior experimental numbers are void and
will be regenerated. This version replaces the heuristic with closed-form theory,
adds a two-tier oracle, error bars, and the restructuring required by advisor feedback.

---

## 1. Problem Setup

`N` equal-radius particles (radius `r`) move in a bounded `d`-dimensional domain
(`d ∈ {2,3}`) of volume `V`. Broad-phase collision detection uses a uniform grid of
cubic cells with edge length `ℓ`. Each particle is binned into its cell; only particles
sharing a cell or a neighboring cell are tested for overlap.

The per-frame cost is a U-shaped function of `ℓ`: small `ℓ` → many cells, cell-traversal
dominates; large `ℓ` → many particles per cell, quadratic intra-cell pair tests dominate.
In **dynamic** scenes (explosions, collapses, periodic oscillation) the cost-optimal `ℓ`
shifts continuously over a run; no single fixed `ℓ` is optimal throughout.

### Notation

| Symbol | Meaning |
|---|---|
| `d` | spatial dimension (2 or 3) |
| `N` | particle count |
| `r` | particle radius (equal for all) |
| `V` | domain volume / area |
| `ℓ` | grid cell edge length |
| `M(ℓ) = V/ℓ^d` | number of cells |
| `n_i` | particle count in cell `i` |
| `S(ℓ) = Σ_i n_i²` | sum of squared cell occupancies |
| `p_i = n_i/N` | occupancy fraction |
| `H₂(ℓ) = −log Σ_i p_i²` | discrete Rényi-2 entropy of the occupancy distribution |
| `h₂` | differential Rényi-2 entropy of the underlying density |
| `D₂` | correlation (Rényi-2) dimension |
| `T(ℓ)` | per-frame wall-clock cost (frame time) |
| `a`, `b`, `c₀` | per-cell cost, per-pair-test cost, ℓ-independent cost |
| `B(ℓ)=a·M(ℓ)`, `Q(ℓ)=b·S(ℓ)` | per-cell / per-pair-test cost components |
| `ℓ*` | cost-optimal cell size |

---

## 2. Cost Model

One frame's grid-dependent work splits into three sources:
- **Insertion** of `N` particles — cost `∝ N`. Since `N` is fixed, this is **ℓ-independent**.
- **Per-cell** work — clear (memset), prefix sum, and cell iteration in the query —
  all `∝ M(ℓ)`.
- **Pair tests** — candidate overlap tests, `∝ S(ℓ) = Σ n_i²`.

Therefore:

$$T(\ell) = c_0 + a\,M(\ell) + b\,S(\ell)$$

`c₀` absorbs all ℓ-independent work; `a` is the per-cell coefficient; `b` is the
per-candidate-pair-test coefficient. For **minimization over ℓ**, only `a·M + b·S`
matter — `c₀` drops out.

**Note on the "two-weight" critique.** The previous draft was criticized for mapping
"four cost sources to two weights." For the purpose of choosing `ℓ*` this is not an
approximation: every per-cell cost (build memset, prefix, query cell-iteration) is
linear in `M`, every pair cost is linear in `S`, and the `N`-term is constant. The
model collapses *exactly* to two ℓ-dependent coefficients. The only modeling
assumptions are: (i) `a, b` are stable, and (ii) cost is linear in `M` and `S`
(no cache/bandwidth nonlinearity). Both are validated empirically (Section 8).

---

## 3. Exact Entropy Identity

With `p_i = n_i/N`, the discrete Rényi-2 entropy of the cell-occupancy distribution is
`H₂(ℓ) = −log Σ_i p_i²`. Then, by definition:

$$S(\ell) = \sum_i n_i^2 = N^2 \sum_i p_i^2 = N^2\, e^{-H_2(\ell)}$$

This is **exact** — the pair-test cost *is* the collision probability of the
discretized distribution. No approximation.

**Cauchy–Schwarz building block.** `Σ n_i² ≥ (Σ n_i)²/M = N²/M`, with equality iff
the distribution is uniform. Hence `H₂(ℓ) ≤ log M(ℓ)`: a uniform distribution
maximizes entropy and minimizes pair-test cost at fixed `ℓ`.

---

## 4. Scaling Law and the Correlation Dimension

To minimize `T(ℓ)` we need how `S` scales with `ℓ`. In a scaling (mesoscopic) regime:

$$S(\ell) \approx S(\ell_0)\,(\ell/\ell_0)^{D_2}$$

The exponent is **not free**: `Σ p_i²(ℓ)` is the Rényi-2 partition function, and

$$D_2 \equiv \lim_{\ell\to 0}\frac{\log \sum_i p_i^2(\ell)}{\log \ell}$$

is the **correlation dimension** (Rényi-2 dimension; Grassberger–Procaccia,
Hentschel–Procaccia). `D₂ = d` for a space-filling density; `D₂ < d` for a
distribution concentrated on a lower-dimensional set (a curve, shell, sheet, filament).

**This is the conceptual core:** the cost of collision detection is governed by the
fractal/correlation dimension of the instantaneous particle cloud. When particles
collapse onto a structure `D₂` drops; when they fill space `D₂ → d`.

**Status:** this power law is the single modeling pillar of the derivation (the
entropy identity in §3 is exact; the cost-model linearity in §2 is the other
assumption). It is validated empirically (Section 8) and the algorithm is
self-diagnosing on it (Section 7).

---

## 5. Main Result — Closed-Form Optimal Cell Size

The ℓ-dependent cost is `T(ℓ) = c₀ + aV·ℓ^{−d} + bK·ℓ^{D₂}` with `K = S_cur/ℓ_cur^{D₂}`.
Setting `dT/dℓ = 0`:

$$d\,aV\,\ell^{-d-1} = D_2\,bK\,\ell^{D_2-1}$$

Solving and expressing via the *measured* current-frame quantities
`B_cur = a·M_cur`, `Q_cur = b·S_cur`:

$$\boxed{\;\ell^* = \ell_{\text{cur}}\left[\frac{d\cdot B_{\text{cur}}}{D_2\cdot Q_{\text{cur}}}\right]^{\frac{1}{d+D_2}}\;}$$

For `D₂ ≥ 1` the cost is strictly convex, so `ℓ*` is the unique minimizer.

**Reading.** Optimal cell size = current cell size scaled by the `(d+D₂)`-th root of the
cost imbalance. If per-cell cost dominates (grid too fine) the bracket exceeds 1 and `ℓ`
grows; if pair-test cost dominates (too coarse) `ℓ` shrinks. It is a multiplicative,
Newton-like step toward balance.

**Corollary — balance principle.** At the optimum:

$$d\cdot B(\ell^*) = D_2\cdot Q(\ell^*)$$

The optimal grid balances the two cost terms, weighted by their scaling exponents.

**Entropy interpretation.** Since `S = N² e^{−H₂}`, higher Rényi-2 entropy (a more
spread-out cloud) shifts `ℓ*`. In the space-filling case `D₂ = d` this reduces to
`ℓ* ∝ exp(h₂/2d)`: the optimal cell size scales exponentially with the differential
Rényi-2 entropy of the distribution.

---

## 6. Convergence

On a static distribution, work in log-space `x = log ℓ`. Define `s(x) = log S(e^x)`,
so the local correlation dimension is `s'(x) = D₂(x)` and its scale-drift is
`s''(x) = dD₂/d log ℓ`. The update is a map `x_cur ↦ x_next`.

**Proposition.** Near the optimum the update map has derivative
`φ'(x*) = − s''(x*) / [D₂*(d + D₂*)]`, so the error contracts each frame by

$$c = \frac{|s''(x^*)|}{D_2^*\,(d + D_2^*)}$$

Consequences:
1. **Exact power law ⟹ one-shot.** `D₂` scale-invariant ⟹ `s'' = 0` ⟹ `c = 0`: one
   application lands exactly on `ℓ*`. (This is how the formula was derived.)
2. **Convergence condition.** `c < 1 ⟺ |dD₂/d log ℓ| < D₂(d+D₂)` (RHS `≥ 3` in 2D,
   `≥ 4` in 3D). Real distributions have a near scale-invariant `D₂` over the relevant
   range ⟹ small `c` ⟹ convergence in 1–3 frames. The sign is negative: convergence
   proceeds by damped alternation (overshoot–correct).
3. **Dynamic tracking.** If `ℓ*` moves by `Δ` per frame, the steady-state tracking lag
   is `L = Δ / |1 − φ'(x*)|` (recurrence `e_{t+1} = φ'(x*)·e_t − Δ`). Since `|φ'| < 1`
   under the convergence condition, `L` lies in `(Δ/2, ∞)`; in the near-one-shot regime
   (`|φ'| ≈ 0`) `L ≈ Δ` — an *irreducible one-frame lag*, because the update can only
   target where `ℓ*` just was. Fast contraction does not remove this floor but prevents
   the lag from compounding: a method that takes many frames to converge, or adapts only
   every `K` frames, instead lags by a multiple of `Δ`. (Numerically verified;
   see `math_kernel_reference.py`.)

`s''` is measurable from the oracle runs' `log S`–`log ℓ` curves, so `c` can be
*predicted* and checked against observed convergence — a theory-validated-in-experiment
loop.

**Mode A note.** With `D₂ = d` assumed, the update converges, but to a fixed point
offset from `ℓ*` by an amount controlled by `|d − D₂*|`; Mode B (measured `D₂`) removes
this bias.

**Assumptions.** Static/slowly-varying distribution; `S` twice-differentiable in
log-log (mesoscopic regime, away from the `ℓ → 2r` breakdown); within the clamp
interval; converged EMA estimates. Estimation noise gives convergence to a
noise-sized neighborhood of `ℓ*`, not a point. (Full Mode-A algebra to be written
during paper formalization.)

---

## 7. The Adaptive Algorithm

Called once per frame, after the collision pass.

```
HER FRAME  (after the collision pass)

State:      ℓ_cur, EMA weights â, b̂, EMA dimension D̂₂, frame counter t
Constants:  d, V, r, K_D2, ε, γ, EMA rates α_w, α_D

1.  // FREE measurements from this frame's grid
    M_cur ← V / ℓ_cur^d
    S_cur ← Σ n_i²                  // computed during grid build — NO extra pass
    t_M   ← time of M-proportional work   // sum of clear+prefix+cell-iteration regions
    t_S   ← time of pair-test region

2.  // update cost coefficients (EMA)
    â ← α_w·(t_M / M_cur) + (1-α_w)·â
    b̂ ← α_w·(t_S / S_cur) + (1-α_w)·b̂

3.  // refresh correlation dimension
    if Mode B and (t mod K_D2 == 0):
        S₋ ← bin at ℓ_cur·γ⁻¹       // 1 extra O(N) pass
        S₊ ← bin at ℓ_cur·γ         // 1 extra O(N) pass
        D_raw ← least-squares slope of log S vs log ℓ
                over { ℓ_cur·γ⁻¹, ℓ_cur, ℓ_cur·γ }
        D̂₂ ← α_D·D_raw + (1-α_D)·D̂₂
    else if Mode A:
        D̂₂ ← d

4.  // closed-form optimal cell size  (Section 5)
    B_cur ← â·M_cur ;  Q_cur ← b̂·S_cur
    ℓ_raw ← ℓ_cur · [ (d·B_cur) / (D̂₂·Q_cur) ]^(1/(d+D̂₂))

5.  // clamp to the scaling regime  (ℓ_min is also a correctness floor)
    ℓ* ← clip(ℓ_raw, ℓ_min, ℓ_max),   ℓ_min = 2r

6.  // hysteresis: resize only on a meaningful change
    if |ln(ℓ*/ℓ_cur)| > ε:
        GridResize(ℓ*) ;  ℓ_cur ← ℓ*

7.  t ← t+1
```

**Mode A** (D₂ = d, *zero* extra passes) — exact for space-filling scenes, biased on
collapsed ones. **Mode B** (probed D₂, amortized `2/K_D2` extra passes per frame) —
robust everywhere. Both crush the previous draft's 7–9 candidate search on overhead.

### Design decisions (with rationale)

- **Three timer regions** (`t_M`, `t_S`, and an `t_N` for insertion, logged but unused
  in the formula). The previous build/query split could not *isolate* `a` and `b`;
  three regions identify them directly. The same data validates cost-model linearity
  (Pillar i) — one instrumentation, two uses.
- **D₂ via 3-point log-log least-squares slope, EMA-smoothed, probed every `K_D2`
  frames** (not every frame). `D₂` varies slowly, so amortizing the probe is free
  accuracy. `γ ≈ 1.3` keeps probe points far enough for a clean slope, close enough
  to stay local. If the 3 points are not collinear in log-log, the power-law
  assumption is locally violated — the algorithm self-diagnoses Pillar ii.
- **`ℓ_min = 2r`** — principled on three counts: the power law breaks there
  (occupancy → {0,1}); it is the Ericson `ℓ = 2r` heuristic re-entering *as a floor*;
  and **it is a correctness floor** — with a 1-ring half-neighbor stencil, `ℓ ≥ 2r`
  is required for completeness (no missed overlapping pairs). `ℓ_max` = loose safety
  bound (e.g. quarter of the minimum domain extent).
- **ℓ-update every frame (K = 1).** The previous `K = 5` existed only to amortize the
  7–9 probes. The closed form is free, so there is no reason to throttle it.
- **Hysteresis `ε` on the log scale** (relative), default `ε = ln(1.1)` (a 10%
  threshold) — avoids resize thrashing on noise.
- **`ℓ` is not EMA-smoothed** — the formula is a direct solve, not an iterative
  search; stability comes from EMA-smoothed inputs (`â, b̂, D̂₂`) plus hysteresis.

### Default parameters (to be ablated, Section 9)

`α_w = 0.7` (fast weight recovery), `α_D ≈ 0.5` (smooth D₂ noise without lagging real
change), `γ = 1.3`, `K_D2 = 5`, `ε = ln(1.1)`, `K = 1`, `ℓ_min = 2r`,
`ℓ_max = min(domain extent)/4`.

---

## 8. Theoretical Pillars to Validate

The derivation rests on exactly two modeling assumptions. The experiments must test
both before any `ℓ*` claim is trusted.

- **Pillar i — cost-model linearity.** `T`'s ℓ-dependent part is linear in `M` and `S`
  with stable `a, b`. *Test:* from the per-frame timer logs, plot `t_M` vs `M` and
  `t_S` vs `S` across all frames/scenarios; check linearity and `R²`; report where
  cache/bandwidth effects bend the line.
- **Pillar ii — power-law scaling.** `S(ℓ) ∝ ℓ^{D₂}` in the scaling regime. *Test:*
  from the oracle runs (which sweep `S(ℓ)` over all candidate `ℓ` each frame), plot
  `log S` vs `log ℓ`; check local straightness; extract `D₂` and its scale-drift
  `s''` (the latter feeds the §6 convergence prediction).

---

## 9. Experiment Design

### 9.1 Scenarios — stressing the D₂ axis

`D₂` is the centerpiece, so the scenario set is chosen so that one group holds
`D₂ ≈ d` and another drives `D₂ < d`. The contrast is the paper's empirical payoff:
Mode A is fine on the first group and degrades on the second; Mode B handles both.

**2D (~8 scenarios):**

| Scenario | Description | Regime stressed |
|---|---|---|
| stable | uniform, random velocities | static, `D₂ ≈ d` — overhead control |
| collapse | particles converge inward | rising density, `D₂ ≈ d` — a/b balance |
| explosion | central cluster bursts radially | transient thin shell → `D₂ < d`, then fills |
| vortex/spiral | particles pile onto spiral arms | filamentary, `D₂ < d` — strong test |
| funnel | sinusoidal attractor | filament formation |
| pulse | periodic radial spring | fast periodic change — tracking lag |
| multicluster | several merging clusters | spatial non-uniformity at `D₂ ≈ d` |
| (8th, optional) stream | particles funneled through a channel | filamentary flow |

**3D (~5 scenarios):** stable, collapse, explosion, pulse, **settling** (gravity-driven
settling into a 2D layer in a 3D domain → clean `D₂ ≈ 2` dimension-drop test).

Each scenario is a **deterministic function of a seed** (initial positions/velocities).

### 9.2 Baselines / methods compared

- **Fixed `ℓ = 2r`** (Ericson) — the universal default.
- **Fixed `ℓ` = Sigurgeirsson analytical formula** — a real published *method* (not
  just a number). *TODO: extract their exact optimal-ℓ expression when writing.*
- **Static oracle** — best single fixed `ℓ` over the whole run (exhaustive over the
  candidate set). Ceiling on any fixed method.
- **Per-frame oracle** — best `ℓ` for *each* frame (exhaustive per frame). Ceiling on
  any *adaptive* method = the true `ℓ_emp` trajectory. The adaptive method is compared
  primarily against this. (New vs the previous draft, which had only the static oracle.)
- **Blind search** — the previous draft's 7–9 candidate probe, now an **ablation**
  ("adaptation without the analytical update rule").
- **Ours — Mode A** and **Ours — Mode B**.
- **(Optional) sweep-and-prune** — a structurally different broad-phase; secondary
  reference for "why a uniform grid at all," framed as a different regime.

### 9.3 Metrics — all anchored to measured frame time

- **Mean ms/frame** and total wall time (primary).
- **Gap to per-frame oracle (%)** — closeness to best achievable FPS.
- **Gap to static oracle** — `> 100%` means the adaptive method beats the best fixed `ℓ`.
- **Adaptation overhead** — separately measured ms spent in the adaptation logic.
  Mode A `≈ 0`, Mode B tiny, blind search large. Substantiates the "near-zero
  overhead" claim directly.
- **ℓ-tracking error** — `|log(ℓ_chosen / ℓ_emp)|` over time.
- **D₂(t) trajectory** — both a validation artifact (shows `D₂` genuinely drops in
  explosion/vortex/settling) and the basis for correlating Mode-A bias with `|d − D₂|`.

### 9.4 Per-frame instrumentation (logged every frame)

`ms/frame`; the three timers `t_M, t_S, t_N`; `ℓ_cur, M_cur, S_cur, D̂₂, ℓ*, B_cur,
Q_cur`, the balance residual `r`. **Oracle runs additionally log** `S(ℓ)` and `ms` at
*every* candidate `ℓ` → yields `ℓ_emp` per frame, the `log S`–`log ℓ` curve, and `s''`.
The oracle runs thus serve three purposes at once: baseline, validation of both
pillars, and the convergence-constant measurement.

### 9.5 Ablations

1. Mode A vs Mode B — the value of measuring `D₂`.
2. Analytical update vs blind search — the value of the closed form (overhead + accuracy).
3. Sensitivity to `K_D2`, `γ`, `ε`, `α_D`.
4. With vs without the clamp — confirms the `ℓ → 2r` breakdown is real.

### 9.6 Setup

`d ∈ {2,3}`. 2D: `N ∈ {50K, 100K}`; 3D: `N = 50K` plus one larger `N` to show scaling.
2D domain `800 × 600`, 3D domain `200³`. 300 frames in 2D, 100 in 3D (long enough for
several pulse periods). Physics: `r = 0.5`, restitution `0.95`, `Δt = 1/60 s`.
**Single-threaded** measured runs (keeps the cost model honest; the multi-core CPU is
used only to run independent configurations in parallel as separate processes).
Fixed RNG seeds; **multiple seeds per scenario** → report mean ± std (the previous
draft had no error bars). Build: `gcc -O2 -march=native`. Hardware: Ryzen 5 5500U,
no GPU (GPU generalization is an explicit, stated limitation, not a goal).

**Run matrix:** {scenario} × {N} × {method} × {seed}, all deterministic. Raw per-frame
logs saved as CSV; all plots regenerated from logs.

---

## 10. Implementation Plan (for Claude Code)

**Language/build.** C, `gcc -O2 -march=native`, single-threaded measured paths.
Config-driven (JSON), deterministic RNG (specify and seed a named generator, e.g.
PCG or xoshiro256++). Git-versioned. Raw logs → CSV; analysis/plots in Python.

**Grid data structure.** Counting-sort binning: (1) compute cell index per particle,
(2) per-cell counts, (3) prefix sum → cell start offsets, (4) scatter into a sorted
particle array. Query: per cell, intra-cell pairs + half-neighbor cells (4 in 2D,
13 in 3D). Narrow phase: circle/sphere overlap, elastic resolution with restitution.

**Correctness requirement.** The broad phase must return exactly the set of truly
overlapping pairs — verify against an `O(N²)` brute-force check on small `N`. Note
`ℓ ≥ 2r` is required for completeness with the 1-ring stencil (see §7).

**Build phases:**
1. Core 2D simulator: particles + uniform grid broad-phase + narrow-phase + integrator,
   fixed `ℓ`. Correctness checks (pair-set vs brute force; energy/momentum sanity).
2. Instrumentation: three timer regions, full per-frame logging, config + seeds.
3. Scenario generators (all 2D scenarios as deterministic seed functions).
4. Baselines: fixed-`ℓ` variants; oracle harness (sweeps all candidate `ℓ`);
   per-frame oracle; blind search; optional sweep-and-prune.
5. The adaptive method: Modes A and B (Section 7).
6. 3D extension (grid, stencil, scenarios incl. settling).
7. Experiment harness: runs the full matrix, collects logs.
8. Analysis: cost-model linearity (Pillar i), power-law check (Pillar ii),
   `D₂(t)` trajectories, predicted vs observed convergence rate, all metrics & plots.

---

## 11. Novelty and Related Work

**Novelty claim.** No prior work derives the optimal uniform-grid cell size in closed
form as an information-theoretic functional of the *actual current* distribution, nor
identifies the correlation dimension as the governing quantity, nor provides a runtime
update rule with the balance principle. The result is also *continuous* (it removes
the previous draft's "discrete candidates" limitation).

**Must engage** (from the previous draft): Ericson (`ℓ=2r`), Green (GPU particles),
Sigurgeirsson et al. (static scaling laws), Apte et al. (a priori search-box sizing),
Teschner et al. and Ihmsen et al. (spatial hashing), Eitz & Gu (hierarchical hashing),
Rényi (entropy). Add Grassberger–Procaccia / Hentschel–Procaccia for the correlation-
dimension definition.

**TODO before/while writing (also addresses advisor point 5):**
- Web-search recent (2021–2025) broad-phase / spatial-data-structure / GPU-particle
  work; add 1–2 recent references and cite them in-text.
- **Verify the novelty claim** against the literature — has anyone connected
  correlation/Rényi dimension to spatial-data-structure performance? If so, reframe.

---

## 12. Open Items and Risks

- **Power law may not hold cleanly** in some scenarios. Either way it is a reported
  result; Mode B's local probe mitigates, and the algorithm self-diagnoses.
- **Cost-model linearity** may bend at very large `M` (cache/bandwidth). If so,
  document as a limitation or add a term.
- **D₂ estimation noise** — mitigated by EMA + probe interval; quantify it.
- **Sigurgeirsson baseline** needs their exact optimal-`ℓ` formula extracted from the
  source paper for a fair comparison.
- **3D cost** is harsher (`M ∝ ℓ^{−3}`); expect the static-overhead story to differ
  from 2D. The closed form handles `d` correctly, but report 3D honestly.
- **Novelty claim** must survive the literature search.

---

## 13. Mapping to Paper Structure

Restructured per advisor feedback. Spec → paper section:

| Paper section | Spec source | Advisor point addressed |
|---|---|---|
| Abstract | §0 core claim | **(1)** one-sentence contribution vs literature |
| Introduction | §0–§1, contributions list | — |
| Related Work — *before* Method | §11 | **(4)** related work precedes the method |
| Method | §2–§7; **add a grid + collision figure** | **(2)** method figure; **(3)** all figures/tables referenced in preceding text |
| Experiments | §8–§9; pillar validation, main comparison, ablations | — |
| Comparison with literature | woven into Related Work + Experiments (a short subsection is acceptable) | **(4)** |
| Limitations | §12 + GPU/CPU scope | — |
| Conclusion | §0 | — |
| References | §11 + recent additions | **(5)** add recent references |

Every figure and table must be introduced by name in the preceding paragraph
("Figure 1 shows an example cell and collision configuration ...") — advisor point (3).
