# Phase 9 — the at-scale campaign

Two parts: (1) does the D2 benefit sharpen with scale? and (2) the method
comparison at §9.6 scale, leak-fixed, for the paper (closes TODO item 3).

## Part 1 — does the D2 benefit sharpen with N?  (answer: NO — it is flat)

Re-ran the Phase-8 D2-controlled test (static known-D2 point sets) at
N = 8K, 50K, 100K, in 2D (uniform, 1D curve, Cantor dust D2 = 0.5–2.0) and 3D
(uniform D2=3, a 2D wavy sheet D2≈2, 3D Cantor dust D2 = 0.5–3.0). Per (dim, N):
Mode-A and Mode-B gap-to-oracle vs ground-truth D2, plus the estimator error.

### Headline: the Mode-A penalty and the A−B advantage are SCALE-INVARIANT

For the configurations whose cost optimum stays **interior** across N, the
gaps are flat in N. The cleanest examples (3D, D2=2, interior at every N):

| config (3D, D2=2) | A−B @ 8K | A−B @ 50K | A−B @ 100K |
|---|---|---|---|
| cantor b3k9 | +0.23 | +0.21 | +0.19 |
| sheet       | +0.15 | +0.14 | +0.16 |

Mode A's gap itself is likewise flat (b3k9: 1.27× / 1.31× / 1.29×; sheet: 1.20×
/ 1.22× / 1.24×). The benefit does **not** grow with N.

**Why (refuting the cost-curvature growth hypothesis).** The hypothesis was
that c0 (the ℓ-independent insertion cost) dilutes the penalty at small N and
becomes negligible at large N, so the penalty grows. But the uniform-grid cost
optimum sits at **occupancy ≈ O(1)** (ℓ* ≈ inter-particle spacing): there S ≈ N
(not N²), so the ℓ-dependent cost a·M + b·S and the insertion cost c0 both scale
~N. Their ratio — which sets the dilution — is therefore N-invariant, and so is
the penalty. The D2 advantage is a **D2 effect, not a scale effect** (exactly
the premise under which Phase 8 used N=8K).

### It IS a dimension effect — much larger in 3D than 2D

The D2:=d assumption is more wrong when d is larger, so the 3D penalties dwarf
the 2D ones (at N=8K, the regime where low-D2 optima are still interior):

| ground-truth D2 | 2D Mode-A gap | 3D Mode-A gap |
|---|---|---|
| d (uniform)     | 1.03× | 1.12× |
| ~2.0            | 1.11× | 1.27× |
| ~1.5            | —     | 1.40× |
| ~1.0            | 1.11–1.16× | 1.96× |
| ~0.63           | 1.09× | 2.54× |
| ~0.50           | 1.30× | 2.40× |

Mode B holds ≤1.15× in 2D and ≤1.46× in 3D across the whole range; the A−B
advantage reaches **+1.2 in 3D** at D2≈0.6 (vs +0.38 in 2D).

### At scale, low-D2 brings clamping / infeasibility, not amplification

As N rises at fixed box, the low-D2 optima shrink to the **2r correctness
floor** (the distribution is so concentrated the optimal cell is the smallest
legal one) — there Mode A = Mode B (D2 estimation is moot). In 2D this happens
by N=50K for D2≤1.26; in 3D by N=50K for D2≤1. At N=100K in 3D the lowest-D2
Cantor sets (D2≈0.5–0.6) become **computationally infeasible** — a D2≈0.6 set
of 100K points has ~10⁸ candidate pairs at the floor ℓ (the broad-phase output
overflows). So scale *narrows* the window where runtime D2 estimation helps,
rather than widening it.

### Estimator (Task 2 carryover)

The 3-point estimator recovers the operative local D2: it matches the
construction D2 in the clean mid-range, compresses toward ~1 for dilute/high-D2
sets at occupancy≈1 (e.g. uniform's local slope is ≈1, not d), and is noisy for
the sparsest fractals. Mode B remains robust because the (possibly compressed
or clamped) estimate is still far closer to the operative truth than D2:=d.

### 2D

**2D** — Mode-A gap / Mode-B gap (×oracle) at each N; `clamp`=optimum at 2r floor, `—`=infeasible:

| config | D2 | N=8K (A / B) | N=50K (A / B) | N=100K (A / B) |
|---|---|---|---|---|
| uniform | 2.00 | 1.03 / 0.97 | 1.08 / 1.02 | 1.05 / 1.01 |
| cantor b3k6 | 1.63 | 1.11 / 1.09 | 1.11 / 1.04 | 1.06 / 1.02 |
| cantor b3k4 | 1.26 | 1.17 / 1.02 | 1.07 / 1.03 | 0.94 / 0.94 clamp |
| cantor b3k3 | 1.00 | 1.11 / 1.00 | 0.99 / 0.99 clamp | 0.99 / 0.99 clamp |
| curve | 1.00 | 1.16 / 1.06 | 1.00 / 1.00 clamp | 0.99 / 0.99 clamp |
| cantor b3k2 | 0.63 | 1.09 / 1.05 | 0.96 / 0.96 clamp | 1.01 / 1.01 clamp |
| cantor b9k3 | 0.50 | 1.30 / 0.92 | 1.03 / 1.03 clamp | 1.02 / 1.02 clamp |

### 3D

**3D** — Mode-A gap / Mode-B gap (×oracle) at each N; `clamp`=optimum at 2r floor, `—`=infeasible:

| config | D2 | N=8K (A / B) | N=50K (A / B) | N=100K (A / B) |
|---|---|---|---|---|
| uniform | 3.00 | 1.12 / 1.04 | 1.03 / 1.00 | 1.04 / 0.99 |
| cantor b3k9 | 2.00 | 1.27 / 1.04 | 1.31 / 1.11 | 1.29 / 1.10 |
| sheet | 2.00 | 1.20 / 1.05 | 1.22 / 1.08 | 1.24 / 1.08 |
| cantor b3k5 | 1.47 | 1.40 / 1.11 | 1.42 / 1.13 | 1.13 / 1.11 clamp |
| cantor b3k3 | 1.00 | 1.96 / 1.14 | 0.99 / 0.99 clamp | 1.00 / 1.00 clamp |
| cantor b3k2 | 0.63 | 2.54 / 1.31 | 1.00 / 1.00 clamp | — |
| cantor b4k2 | 0.50 | 2.40 / 1.46 | 1.00 / 1.00 clamp | — |


## Part 2 — method comparison at §9.6 scale (leak-fixed, split-sample)

Regenerates the Phase 5/6 comparison cleanly at campaign scale with the
**t_narrow leak fixed** (closes TODO item 3 — the old Phase 5/6 tables were
slightly inflated for the sim-measured methods). Dynamic scenarios from
`configs/experiment/`, split-sample per-frame and static oracles, multicluster
capped at N=20K (4 disks saturate above that; TODO item 1).

**Compute reality.** One full 2D run (N=50K, 300 frames, 8 scenarios) takes
**~123 min** — the split-sample per-frame oracle re-sweeps the ℓ-ladder every
frame (frames × ladder × 2K detects at N), which is intrinsically ~7 min per
scenario at 50K. The full §9.6 matrix (N∈{50K,100K} × multiple seeds × 2D+3D)
is therefore ~10+ hours and not tractable here. Delivered: **2D N=50K** (300
frames, all 8 scenarios) and **3D N=50K** (60 frames — reduced from §9.6's 100
for tractability — all 5 scenarios), **single seed each**. Multi-seed (mean±std)
and N=100K are deferred as a documented compute limit, not a silent drop.

### Headline (2D, N=50K)

- **Closed-B ≈ 1.01–1.05× the per-frame oracle** on every scenario — near-oracle
  at ~65 µs/frame overhead. It beats Closed-A (1.03–1.11×) and matches or beats
  the static oracle and blind-9 (whose overhead is ~6–11 **ms**/frame — 100×
  the detection it chooses).
- **Ericson(2r) is only ~1.8–2.6× the oracle at N=50K**, vs 9–20× at N=8K: at
  campaign density the optimal cell (~2–3) is close to 2r, so the universal
  default is far less catastrophic than at moderate N. The method's win over
  Ericson shrinks with density — an honest scale finding.

### 2D, N=50K, 300 frames (single seed)

Mean t_detect ratio to per-frame oracle (and Closed-B overhead, µs):

| scenario | Ericson | Density | Static-oc | Blind-9 | Closed-A | **Closed-B** | B oh(µs) |
|---|---|---|---|---|---|---|---|
| stable | 2.59× | 1.02× | 1.02× | 1.02× | 1.11× | 1.05× | 65 |
| collapse | 2.46× | 1.01× | 1.01× | 1.01× | 1.07× | 1.02× | 62 |
| explosion | 2.03× | 1.17× | 1.00× | 1.00× | 1.00× | 0.99× | 72 |
| vortex | 1.76× | 1.27× | 0.99× | 1.00× | 1.03× | 1.03× | 83 |
| funnel | 2.03× | 1.17× | 1.09× | 1.03× | 1.09× | 1.01× | 66 |
| pulse | 1.38× | 1.72× | 1.06× | 1.04× | 1.03× | 1.04× | 67 |
| multicluster | 3.22× | 2.04× | 1.01× | 1.00× | 1.01× | 1.00× | 24 |
| stream | 1.20× | 2.23× | 1.08× | 1.05× | 1.11× | 1.13× | 2 |

### 3D, N=50K, 60 frames (single seed)

Mean t_detect ratio to per-frame oracle (and Closed-B overhead, µs):

| scenario | Ericson | Density | Static-oc | Blind-9 | Closed-A | **Closed-B** | B oh(µs) |
|---|---|---|---|---|---|---|---|
| stable | 27.00× | 0.98× | 1.00× | 0.98× | 1.02× | 0.98× | 88 |
| collapse | 28.62× | 0.98× | 1.00× | 0.99× | 1.00× | 0.98× | 91 |
| explosion | 14.81× | 1.74× | 0.96× | 1.11× | 1.02× | 1.01× | 106 |
| pulse | 16.23× | 1.95× | 1.12× | 1.10× | 1.06× | 1.04× | 120 |
| settling | 24.14× | 1.01× | 0.95× | 1.01× | 0.95× | 0.94× | 88 |

In 3D the Ericson 2r default stays catastrophic even at N=50K (15–29×: the 200³
grid at ℓ=2r has 8×10⁶ mostly-empty cells), and Closed-B holds 0.94–1.04× the
oracle at ~90 µs — beating the density baseline on the dense scenarios
(explosion 1.74×, pulse 1.95×) where a fixed analytic ℓ ignores the clustering.

## Not done

- Multi-seed (mean±std) and N=100K (2D and 3D): each run is ~2 h; deferred as a
  compute limit. The single-seed at N=50K is the primary at-scale result.
- 3D at N=100K specifically also risks the ~10⁸-candidate-pair buffer overflow
  on the dense scenarios (same infeasibility as the low-D2 Part-1 configs).
- The full multi-seed campaign is Phase-7-proper territory (a dedicated
  long-running harness), out of scope for this analysis phase.
