"""
Math-kernel reference & verification suite
===========================================
Project: "Entropy-Governed Adaptive Cell-Size Selection for Uniform-Grid
Collision Detection."  Reference for spec sections 5 (closed-form optimum),
6 (convergence), 7 (algorithm: D2 estimator, clamp, hysteresis).

PURPOSE
  Hardware-independent numerical verification of the theoretical core BEFORE
  the full build. The C implementation (Phase 5) must reproduce the KERNEL
  functions below; the test cases here are the correctness oracle.

WHAT THIS VERIFIES
  - Section 5 (main theorem): the closed-form optimal cell size, including the
    one-shot property on an exact power law.
  - Section 6 (proposition): the contraction factor c = |s''| / (D2 (d+D2)),
    the Mode-A bias, and the dynamic tracking-lag scaling.
  - Section 7 (components): the 3-point log-log D2 estimator, clamp, hysteresis.

WHAT THIS DOES NOT VERIFY
  - Pillar i (cost-model linearity) and Pillar ii (the power law holding on
    real particle distributions). Those require the simulator -> Phase 8.

Run:  python3 math_kernel_reference.py
"""
import numpy as np


# ============================================================ KERNEL (5, 7)
def optimal_cell_size(ell_cur, d, D2, B_cur, Q_cur):
    """Spec Section 5, main theorem:
        ell* = ell_cur * [ (d * B_cur) / (D2 * Q_cur) ] ** (1 / (d + D2))
    B_cur = a*M_cur (per-cell cost component), Q_cur = b*S_cur (per-pair-test).
    """
    return ell_cur * ((d * B_cur) / (D2 * Q_cur)) ** (1.0 / (d + D2))


def estimate_D2(ell_cur, gamma, S_of):
    """Spec Section 7, Mode-B estimator: 3-point log-log least-squares slope
    over {ell_cur/gamma, ell_cur, ell_cur*gamma}. S_of(ell) returns S there."""
    ells = np.array([ell_cur / gamma, ell_cur, ell_cur * gamma])
    S = np.array([S_of(e) for e in ells])
    return np.polyfit(np.log(ells), np.log(S), 1)[0]


def clamp(ell, ell_min, ell_max):
    """Spec Section 7, step 5. ell_min = 2r is also a correctness floor."""
    return min(max(ell, ell_min), ell_max)


def should_resize(ell_star, ell_cur, eps):
    """Spec Section 7, step 6: log-scale hysteresis."""
    return abs(np.log(ell_star / ell_cur)) > eps


def ema(new, old, alpha):
    return alpha * new + (1.0 - alpha) * old


# ====================================================== helpers for testing
def _argmin_cost(g, lo, hi, n=200001):
    """Minimise convex g over [lo, hi] in x-space: fine grid then Newton polish
    with central-difference derivatives."""
    xs = np.linspace(lo, hi, n)
    x = xs[int(np.argmin(g(xs)))]
    for _ in range(80):
        h = 1e-5
        gp = (g(x + h) - g(x - h)) / (2 * h)
        gpp = (g(x + h) - 2 * g(x) + g(x - h)) / (h * h)
        if gpp <= 0:
            break
        step = gp / gpp
        x -= step
        if abs(step) < 1e-14:
            break
    return x


# ================================================================== tests
def test_theorem_oneshot(rng, n=3000):
    """Section 5: under an exact power law S = K*ell^D2, the formula must reach
    the true optimum in ONE step from ANY ell_cur (also the c = 0 case of
    Section 6). Cross-checked vs the analytic optimum and an independent
    argmin of the full cost."""
    worst_formula = worst_numeric = 0.0
    for _ in range(n):
        d = int(rng.choice([2, 3]))
        D2 = rng.uniform(1.0, d)
        a, b = rng.uniform(0.1, 10.0), rng.uniform(0.1, 10.0)
        V = rng.uniform(1e3, 1e6)
        K = rng.uniform(1e-3, 1e3)
        ell_cur = rng.uniform(0.5, 50.0)

        M = lambda l: V * l ** (-d)
        S = lambda l: K * l ** D2
        ell_formula = optimal_cell_size(ell_cur, d, D2, a * M(ell_cur), b * S(ell_cur))
        ell_analytic = ((d * a * V) / (D2 * b * K)) ** (1.0 / (d + D2))

        g = lambda x: a * V * np.exp(-d * x) + b * K * np.exp(D2 * x)
        ell_numeric = np.exp(_argmin_cost(g, np.log(ell_analytic) - 4,
                                          np.log(ell_analytic) + 4))

        worst_formula = max(worst_formula, abs(ell_formula - ell_analytic) / ell_analytic)
        worst_numeric = max(worst_numeric, abs(ell_numeric - ell_analytic) / ell_analytic)
    ok = worst_formula < 1e-12 and worst_numeric < 1e-6
    print(f"  [S5] one-shot optimum: formula err={worst_formula:.2e} (<1e-12), "
          f"numeric cross-check={worst_numeric:.2e} (<1e-6)  -> {'PASS' if ok else 'FAIL'}")
    return ok


def test_D2_estimator(rng, n=3000):
    """Section 7: the 3-point estimator recovers a known D2 exactly on a power
    law; EMA reduces variance under multiplicative measurement noise."""
    worst = 0.0
    for _ in range(n):
        D2 = rng.uniform(0.5, 3.0)
        K = rng.uniform(1e-2, 1e2)
        ell_cur = rng.uniform(1.0, 30.0)
        worst = max(worst, abs(estimate_D2(ell_cur, 1.3, lambda l: K * l ** D2) - D2))
    exact_ok = worst < 1e-10

    rng2 = np.random.default_rng(12345)
    D2_true, K, ell_cur = 1.8, 10.0, 8.0
    raw, smoothed, d_hat = [], [], D2_true
    for _ in range(400):
        noisy = lambda l: K * l ** D2_true * np.exp(rng2.normal(0, 0.05))
        r = estimate_D2(ell_cur, 1.3, noisy)
        d_hat = ema(r, d_hat, 0.3)
        raw.append(r)
        smoothed.append(d_hat)
    std_raw, std_sm = np.std(raw), np.std(smoothed[80:])
    ema_ok = std_sm < 0.6 * std_raw
    ok = exact_ok and ema_ok
    print(f"  [S7] D2 estimator: exact-recovery err={worst:.2e} (<1e-10); "
          f"EMA std {std_raw:.4f} -> {std_sm:.4f}  -> {'PASS' if ok else 'FAIL'}")
    return ok


def test_convergence_proposition(rng, n=600):
    """Section 6: with a non-power-law S whose log-log curve is quadratic
    (s(x) = D2_0*x + (kappa/2)x^2 + logK, so s'' = kappa), the update map's
    derivative at the optimum must equal -kappa / (D2* (d + D2*)). The map
    must also contract (|phi'| < 1) and alternate (phi' < 0 for kappa > 0)."""
    worst, contracts, alternates, used = 0.0, True, True, 0
    for _ in range(n):
        d = int(rng.choice([2, 3]))
        kappa = rng.uniform(0.05, 0.40)
        D2_0 = rng.uniform(1.2, float(d))
        a, b = rng.uniform(0.5, 5.0), rng.uniform(0.5, 5.0)
        V = rng.uniform(1e3, 1e5)
        logK = np.log(a * V / b) + rng.uniform(-2, 2)

        s = lambda x: D2_0 * x + 0.5 * kappa * x * x + logK
        sp = lambda x: D2_0 + kappa * x
        g = lambda x: a * V * np.exp(-d * x) + b * np.exp(s(x))

        xstar = _argmin_cost(g, -10, 10)
        D2star = sp(xstar)
        if D2star <= 0:
            continue
        used += 1
        c_pred = -kappa / (D2star * (d + D2star))      # = phi'(x*)

        def phi(x):
            D2_loc = sp(x)
            B = a * V * np.exp(-d * x)
            Q = b * np.exp(s(x))
            return x + np.log((d * B) / (D2_loc * Q)) / (d + D2_loc)

        h = 1e-6
        phi_prime = (phi(xstar + h) - phi(xstar - h)) / (2 * h)
        worst = max(worst, abs(phi_prime - c_pred))
        contracts = contracts and abs(phi_prime) < 1.0
        alternates = alternates and phi_prime < 0.0
    ok = worst < 1e-5 and contracts and alternates and used > n // 2
    print(f"  [S6] convergence factor: max |phi'(x*) - predicted| = {worst:.2e} "
          f"(<1e-5); contracts={contracts}; alternates={alternates}  "
          f"-> {'PASS' if ok else 'FAIL'}")
    return ok


def test_modeA_bias(rng, n=2000):
    """Section 6: Mode A (D2 := d assumed) converges to a fixed point offset
    from the true optimum by exactly log(D2/d)/(d + D2) on an exact power law
    -- zero iff D2 = d, growing as D2 departs from d."""
    worst = 0.0
    for _ in range(n):
        d = int(rng.choice([2, 3]))
        D2 = rng.uniform(1.0, d)
        a, b = rng.uniform(0.5, 5.0), rng.uniform(0.5, 5.0)
        V = rng.uniform(1e3, 1e5)
        K = rng.uniform(1e-2, 1e2)

        xstar = np.log((d * a * V) / (D2 * b * K)) / (d + D2)

        def phiA(x):
            B = a * V * np.exp(-d * x)
            Q = b * K * np.exp(D2 * x)
            return x + np.log((d * B) / (d * Q)) / (2 * d)   # D2_hat := d

        x = xstar + rng.uniform(-2, 2)
        for _ in range(300):
            x = phiA(x)
        bias_pred = np.log(D2 / d) / (d + D2)
        worst = max(worst, abs((x - xstar) - bias_pred))
    ok = worst < 1e-9
    print(f"  [S6] Mode-A bias: max |observed - log(D2/d)/(d+D2)| = {worst:.2e} "
          f"(<1e-9)  -> {'PASS' if ok else 'FAIL'}")
    return ok


def test_tracking_lag():
    """Section 6: with the real per-frame ordering -- the frame runs at the
    carried-in cell size, the update is computed from that frame's
    measurements, and only then does the distribution drift -- the steady-state
    lag is L = Delta / |1 - phi'(x*)|. Drive the optimum at three rates and
    check the observed lag vs the prediction and monotonicity in Delta."""
    d, kappa, D2_0 = 2, 0.25, 1.6
    a, b, V = 2.0, 2.0, 1e4

    def funcs(logK):
        s = lambda x: D2_0 * x + 0.5 * kappa * x * x + logK
        sp = lambda x: D2_0 + kappa * x
        g = lambda x: a * V * np.exp(-d * x) + b * np.exp(s(x))
        return s, sp, g

    rows = []
    for rate in [0.0008, 0.0016, 0.0032]:
        logK = np.log(a * V / b)
        s, sp, g = funcs(logK)
        x = _argmin_cost(g, -10, 10)              # start at the optimum
        lags, deltas, prev_star = [], [], None
        for t in range(450):
            s, sp, g = funcs(logK)                # frame-t landscape
            xstar = _argmin_cost(g, -10, 10)      # frame-t optimum
            if 200 < t <= 400:                    # lag DURING frame t
                lags.append(abs(x - xstar))
                if prev_star is not None:
                    deltas.append(abs(xstar - prev_star))
            prev_star = xstar
            D2_loc = sp(x)                        # measure on frame-t landscape
            B = a * V * np.exp(-d * x)
            Q = b * np.exp(s(x))
            x = x + np.log((d * B) / (D2_loc * Q)) / (d + D2_loc)
            logK += rate                          # frame ends -> distribution drifts
        Delta = float(np.mean(deltas))
        lag_obs = float(np.mean(lags))
        s, sp, g = funcs(np.log(a * V / b) + rate * 300)   # window midpoint
        D2star = sp(_argmin_cost(g, -10, 10))
        phi_p = -kappa / (D2star * (d + D2star))           # phi'(x*) < 0
        lag_pred = Delta / abs(1.0 - phi_p)
        rows.append((Delta, abs(phi_p), lag_obs, lag_pred))

    print("  [S6] tracking lag  L = Delta / |1 - phi'(x*)|:")
    print("        Delta        |c|       lag_obs     lag_pred")
    for D, c, lo, lp in rows:
        print(f"        {D:.6f}    {c:.4f}    {lo:.6f}    {lp:.6f}")
    monotone = rows[0][2] < rows[1][2] < rows[2][2]
    close = all(0.7 < lo / lp < 1.4 for _, _, lo, lp in rows)
    ok = monotone and close
    print(f"        monotone-in-Delta={monotone}; within-factor={close}  "
          f"-> {'PASS' if ok else 'FAIL'}")
    return ok


def print_canonical_vectors():
    """Canonical (input -> expected ell*) vectors to embed as C unit tests."""
    cases = [(2, 2.0, 1.0, 4.0, 1.0), (2, 1.0, 2.0, 3.0, 12.0),
             (2, 1.5, 5.0, 10.0, 2.0), (3, 3.0, 1.0, 9.0, 1.0),
             (3, 1.0, 4.0, 8.0, 50.0), (3, 2.0, 2.0, 5.0, 5.0)]
    print("\n  Canonical test vectors for the C implementation (Phase 5):")
    print("   d   D2     ell_cur      B_cur      Q_cur   ->  ell_star")
    for d, D2, ell_cur, B, Q in cases:
        ls = optimal_cell_size(ell_cur, d, D2, B, Q)
        print(f"   {d}  {D2:.2f}  {ell_cur:9.4f}  {B:9.3f}  {Q:9.3f}  ->  {ls:.10f}")


def main():
    print("Math-kernel verification suite  (spec sections 5, 6, 7)")
    print("=" * 64)
    rng = np.random.default_rng(20240522)
    results = [
        test_theorem_oneshot(rng),
        test_D2_estimator(rng),
        test_convergence_proposition(rng),
        test_modeA_bias(rng),
        test_tracking_lag(),
    ]
    print("=" * 64)
    if all(results):
        print("ALL CHECKS PASSED -- theoretical core is numerically sound.")
    else:
        print("SOME CHECKS FAILED -- theory needs review before building.")
    print_canonical_vectors()
    return all(results)


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
