#ifndef ADAPTER_H
#define ADAPTER_H

#include <stdint.h>

/*
 * Phase 5: the entropy-governed adaptive cell-size controller (spec §7).
 *
 * Implements the closed-form update
 *
 *     ell* = ell_cur * [ (d * B_cur) / (D2 * Q_cur) ] ^ (1 / (d + D2))
 *
 * (§5 main theorem), the 3-point log-log D2 estimator, EMA smoothing
 * of a, b, D2, log-scale hysteresis, and the clamps [ell_min, ell_max]
 * for ell and [d2_floor, d] for D2.
 *
 *   Mode A: D2 := d (zero extra passes -- exact for D2 = d, biased on
 *           collapsed distributions; the bias is log(D2/d)/(d+D2)).
 *   Mode B: D2 estimated from a 3-point log-log probe at
 *           { ell_cur/gamma, ell_cur, ell_cur*gamma }, every K_D2 frames,
 *           EMA-smoothed. The probe is COUNT-ONLY (grid_count_S_at) so
 *           the cost is two O(N) passes amortised over K_D2 frames.
 *
 * Determinism: the adapter has no internal RNG. Given the same
 * sequence of (t_M, t_S, S_cur, M_cur, snapshot) the controller's
 * floating-point output ell_next is bit-identical across runs *modulo*
 * the timer inputs (which are wall-clock; in pure tests we drive the
 * adapter with synthetic timers and the output is exactly reproducible).
 *
 * The math kernels (adapter_optimal_ell, adapter_estimate_D2, clamp,
 * hysteresis) reproduce reference/math_kernel_reference.py to within
 * 1e-12 relative on the canonical test vectors (see tests/test_phase5.c).
 */

typedef enum {
    ADAPTER_MODE_OFF    = 0,    /* controller disabled; ell fixed       */
    ADAPTER_MODE_A      = 1,    /* D2 := d (zero overhead)              */
    ADAPTER_MODE_B      = 2,    /* D2 from 3-point probe                */
    ADAPTER_MODE_BLIND  = 3,    /* old-style 7-9 candidate search       */
} adapter_mode_t;

typedef struct {
    /* configuration */
    adapter_mode_t mode;
    int      d;                 /* spatial dimension (2 here)           */
    double   V;                 /* domain volume / area                 */
    double   domain_w;
    double   domain_h;
    double   ell_min;           /* = 2r                                 */
    double   ell_max;
    double   gamma;             /* D2-probe stencil ratio, e.g. 1.3     */
    int      K_D2;              /* probe every K_D2 frames              */
    double   alpha_w;           /* EMA rate for a, b                    */
    double   alpha_D;           /* EMA rate for D2                      */
    double   eps;               /* log-scale hysteresis threshold       */
    double   d2_floor;          /* lower clamp on D2 (away from 0)      */
    double   d2_ceil;           /* upper clamp on D2 (== d)             */

    /* persistent state */
    double   ell_cur;
    double   a_hat, b_hat;
    double   D2_hat;
    int      have_weights;
    int      have_D2;
    int64_t  frame;

    /* diagnostics from the most recent adapter_step() */
    double   last_ell_raw;      /* before clamping/hysteresis           */
    double   last_ell_next;     /* after clamping; written even if      */
                                /* hysteresis suppressed the resize     */
    int      last_resized;      /* 1 iff |log(ell_star/ell_cur)| > eps  */
    double   last_D2_raw;       /* most recent probe slope (Mode B);    */
                                /* 0 if no probe this frame             */
    int      last_did_probe;    /* 1 iff D2 was probed this frame       */
    double   last_overhead_s;   /* adapter step wall time (s)           */

    /* blind-search bookkeeping */
    int      blind_period;
    double   blind_gamma;
    double   blind_radius;      /* needed to evaluate candidates        */
    int      blind_K;           /* number of candidates (~9)            */
} adapter_t;


/* ---- math kernels: reproduce math_kernel_reference.py exactly --- */

/* ell* = ell_cur * [ (d*B) / (D2*Q) ] ^ (1/(d+D2)). */
double adapter_optimal_ell(double ell_cur, int d, double D2,
                           double B_cur, double Q_cur);

/* 3-point log-log least-squares slope of log S vs log ell over the stencil
 * { ell_lo, ell_mid, ell_hi }. All S must be > 0. Matches np.polyfit. */
double adapter_estimate_D2(double ell_lo, double ell_mid, double ell_hi,
                           double S_lo, double S_mid, double S_hi);

/* clip(x, lo, hi). */
static inline double adapter_clamp(double x, double lo, double hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

/* True iff |log(ell_star / ell_cur)| > eps. */
int adapter_should_resize(double ell_star, double ell_cur, double eps);

/* EMA: alpha * new + (1-alpha) * old. */
static inline double adapter_ema(double new_v, double old_v, double alpha) {
    return alpha * new_v + (1.0 - alpha) * old_v;
}


/* ---- controller --- */

void adapter_init_default(adapter_t *a, adapter_mode_t mode,
                          int d, double domain_w, double domain_h,
                          double ell_start, double radius);

/* One control step.
 *   t_M, t_S:   measured wall-clock (s) for the M-proportional work and
 *               the pair-test work (broad + narrow) at the CURRENT ell.
 *   S_cur, M_cur: measured S, M at the current ell.
 *   n, px, py:  current particle state (Mode B only -- may be 0/NULL in
 *               Mode A or when probes are skipped). */
double adapter_step(adapter_t *a,
                    double t_M, double t_S,
                    long S_cur, int M_cur,
                    int n, const double *px, const double *py);

/* Blind-search baseline. Evaluates ~blind_K candidate ells around ell_cur
 * (geometric, ratio blind_gamma) every blind_period frames and picks the
 * argmin of full broad+narrow t_detect. Returns the chosen next-frame ell.
 * Internally rebuilds a fresh grid per candidate -- the expensive path the
 * closed-form method is meant to obsolete. */
double adapter_blind_step(adapter_t *a,
                          int n, const double *px, const double *py);

#endif
