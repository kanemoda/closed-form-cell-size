#define _POSIX_C_SOURCE 199309L
#include "adapter.h"
#include "grid.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static inline double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* ============================================================ math kernels */

double adapter_optimal_ell(double ell_cur, int d, double D2,
                           double B_cur, double Q_cur) {
    /* Caller must clamp D2 > 0 (formula has D2 in the denominator); B, Q,
     * and ell_cur must all be > 0. Returns ell_cur unchanged on any
     * degenerate input so a one-off bad frame doesn't NaN the trajectory. */
    if (D2 <= 0.0 || Q_cur <= 0.0 || B_cur <= 0.0 || ell_cur <= 0.0 || d <= 0) {
        return ell_cur;
    }
    double ratio = (((double)d) * B_cur) / (D2 * Q_cur);
    double expt  = 1.0 / ((double)d + D2);
    return ell_cur * pow(ratio, expt);
}

double adapter_estimate_D2(double ell_lo, double ell_mid, double ell_hi,
                           double S_lo, double S_mid, double S_hi) {
    if (S_lo <= 0.0 || S_mid <= 0.0 || S_hi <= 0.0) return 0.0;
    /* 3-point least-squares slope of y = log S vs x = log ell. */
    double x0 = log(ell_lo), x1 = log(ell_mid), x2 = log(ell_hi);
    double y0 = log(S_lo),   y1 = log(S_mid),   y2 = log(S_hi);
    double xbar = (x0 + x1 + x2) / 3.0;
    double ybar = (y0 + y1 + y2) / 3.0;
    double dx0 = x0 - xbar, dx1 = x1 - xbar, dx2 = x2 - xbar;
    double dy0 = y0 - ybar, dy1 = y1 - ybar, dy2 = y2 - ybar;
    double num = dx0 * dy0 + dx1 * dy1 + dx2 * dy2;
    double den = dx0 * dx0 + dx1 * dx1 + dx2 * dx2;
    if (den <= 0.0) return 0.0;
    return num / den;
}

int adapter_should_resize(double ell_star, double ell_cur, double eps) {
    if (ell_cur <= 0.0 || ell_star <= 0.0) return 0;
    double r = log(ell_star / ell_cur);
    return fabs(r) > eps;
}

/* ==================================================== cost-model regression */

cost_fit_t adapter_fit_cost_model(const double *M, const double *S,
                                  const double *t, int n) {
    cost_fit_t f;
    memset(&f, 0, sizeof(f));
    f.n = n;
    if (n < 3 || !M || !S || !t) return f;   /* 3 params: need >= 3 points */

    /* Predictor / response means (center to drop the intercept). */
    double Mbar = 0.0, Sbar = 0.0, tbar = 0.0;
    for (int i = 0; i < n; i++) { Mbar += M[i]; Sbar += S[i]; tbar += t[i]; }
    Mbar /= n; Sbar /= n; tbar /= n;

    /* Centered cross-products. */
    double Smm = 0, Sss = 0, Sms = 0, Smt = 0, Sst = 0, Stt = 0;
    for (int i = 0; i < n; i++) {
        double dm = M[i] - Mbar, ds = S[i] - Sbar, dt = t[i] - tbar;
        Smm += dm * dm; Sss += ds * ds; Sms += dm * ds;
        Smt += dm * dt; Sst += ds * dt; Stt += dt * dt;
    }

    /* Standardise predictors -> the 2x2 normal matrix is [n, n*r; n*r, n]
     * (r = corr(M,S)); det = n^2 (1 - r^2). Keeps magnitudes O(1) despite
     * M ~ 1e4 cells and S ~ 1e5. Degenerate if either predictor is constant
     * (no spread) or M, S are perfectly collinear. */
    double sdM = sqrt(Smm / (double)n);
    double sdS = sqrt(Sss / (double)n);
    if (sdM <= 0.0 || sdS <= 0.0) return f;   /* a predictor has no spread */

    double Suu = Smm / (sdM * sdM);           /* = n */
    double Sww = Sss / (sdS * sdS);           /* = n */
    double Suw = Sms / (sdM * sdS);
    double Sut = Smt / sdM;
    double Swt = Sst / sdS;
    double det = Suu * Sww - Suw * Suw;
    if (det <= 0.0) return f;                 /* collinear M, S */

    double au = (Sut * Sww - Swt * Suw) / det;   /* standardised coeffs */
    double bw = (Suu * Swt - Suw * Sut) / det;
    double a  = au / sdM;                        /* back to raw units */
    double b  = bw / sdS;
    double c0 = tbar - a * Mbar - b * Sbar;
    f.a = a; f.b = b; f.c0 = c0;

    /* Residuals + R^2 against the raw (un-standardised) data. */
    double ss_res = 0, sum_abs = 0, max_abs = 0, sum_rel = 0;
    for (int i = 0; i < n; i++) {
        double pred = c0 + a * M[i] + b * S[i];
        double r = t[i] - pred, ar = fabs(r);
        ss_res += r * r;
        sum_abs += ar;
        if (ar > max_abs) max_abs = ar;
        if (t[i] > 0.0) sum_rel += ar / t[i];
    }
    f.r2             = (Stt > 0.0) ? (1.0 - ss_res / Stt) : 0.0;
    f.rmse           = sqrt(ss_res / (double)n);
    f.mean_abs_resid = sum_abs / (double)n;
    f.max_abs_resid  = max_abs;
    f.mean_rel_resid = sum_rel / (double)n;
    f.ok             = 1;
    return f;
}

void adapter_set_calibrated_weights(adapter_t *a, double a_coef, double b_coef) {
    a->a_hat        = a_coef;
    a->b_hat        = b_coef;
    a->have_weights = 1;
    a->calibrated   = 1;
}

/* ============================================================ controller */

void adapter_init_default(adapter_t *a, adapter_mode_t mode,
                          int d, double domain_w, double domain_h,
                          double ell_start, double radius) {
    memset(a, 0, sizeof(*a));
    a->mode     = mode;
    a->d        = d;
    a->domain_w = domain_w;
    a->domain_h = domain_h;
    /* V is informational only -- the formula uses the measured M_cur, not V.
     * The 3D caller sets a->domain_d after init (needed by the Mode-B probe). */
    a->V        = domain_w * domain_h;
    a->ell_cur  = ell_start;
    a->ell_min  = 2.0 * radius;
    /* ell_max heuristic: a quarter of the minimum domain extent (spec §7). */
    a->ell_max  = 0.25 * fmin(domain_w, domain_h);
    a->gamma    = 1.3;
    a->K_D2     = 5;
    a->alpha_w  = 0.7;
    a->alpha_D  = 0.5;
    a->eps      = log(1.1);
    a->d2_floor = 0.2;
    a->d2_ceil  = (double)d;

    a->D2_hat        = (double)d;
    a->have_weights  = 0;
    a->have_D2       = (mode == ADAPTER_MODE_A);
    a->frame         = 0;

    a->blind_period  = 5;
    a->blind_gamma   = 1.25;
    a->blind_radius  = radius;
    a->blind_K       = 9;
}

double adapter_step(adapter_t *a,
                    double t_M, double t_S,
                    long S_cur, int M_cur,
                    int n, const double *px, const double *py, const double *pz) {
    double t_start = now_s();
    a->last_did_probe = 0;
    a->last_D2_raw    = 0.0;

    if (a->mode == ADAPTER_MODE_OFF) {
        a->last_ell_raw    = a->ell_cur;
        a->last_ell_next   = a->ell_cur;
        a->last_resized    = 0;
        a->last_overhead_s = now_s() - t_start;
        a->frame++;
        return a->ell_cur;
    }

    /* (2) update cost coefficients via EMA -- UNLESS calibrated. When the
     * coefficients come from the one-time regression (adapter_set_calibrated_
     * weights), a, b are fixed hardware constants and the per-frame sub-timers
     * are NOT used to update them: the t_M/t_S split mis-attributes the
     * broad-phase per-cell loop overhead into t_S and inflates b. The
     * sub-timers still arrive here (and are logged by callers) for the
     * paper's Pillar-i analysis -- they just no longer drive a, b. */
    if (!a->calibrated && M_cur > 0 && S_cur > 0) {
        double a_new = t_M / (double)M_cur;
        double b_new = t_S / (double)S_cur;
        if (!a->have_weights) {
            a->a_hat = a_new;
            a->b_hat = b_new;
            a->have_weights = 1;
        } else {
            a->a_hat = adapter_ema(a_new, a->a_hat, a->alpha_w);
            a->b_hat = adapter_ema(b_new, a->b_hat, a->alpha_w);
        }
    }

    /* (3) refresh D2. */
    if (a->mode == ADAPTER_MODE_A) {
        a->D2_hat  = (double)a->d;
        a->have_D2 = 1;
    } else if (a->mode == ADAPTER_MODE_B) {
        if ((a->frame % (int64_t)a->K_D2) == 0 && n > 0 && px && py) {
            double ell_lo = a->ell_cur / a->gamma;
            double ell_hi = a->ell_cur * a->gamma;
            /* Only probe if the stencil is interior to [ell_min, ell_max]. */
            if (ell_lo >= a->ell_min && ell_hi <= a->ell_max
                && ell_lo < a->ell_cur && a->ell_cur < ell_hi) {
                long S_lo, S_hi;
                if (a->d == 3) {
                    S_lo = grid_count_S_at3d(a->domain_w, a->domain_h, a->domain_d,
                                             n, px, py, pz, ell_lo);
                    S_hi = grid_count_S_at3d(a->domain_w, a->domain_h, a->domain_d,
                                             n, px, py, pz, ell_hi);
                } else {
                    S_lo = grid_count_S_at(a->domain_w, a->domain_h, n, px, py, ell_lo);
                    S_hi = grid_count_S_at(a->domain_w, a->domain_h, n, px, py, ell_hi);
                }
                if (S_lo > 0 && S_hi > 0 && S_cur > 0) {
                    double D_raw = adapter_estimate_D2(ell_lo, a->ell_cur, ell_hi,
                                                       (double)S_lo, (double)S_cur,
                                                       (double)S_hi);
                    a->last_D2_raw    = D_raw;
                    a->last_did_probe = 1;
                    if (!a->have_D2) {
                        a->D2_hat = D_raw;
                        a->have_D2 = 1;
                    } else {
                        a->D2_hat = adapter_ema(D_raw, a->D2_hat, a->alpha_D);
                    }
                }
            }
        }
    }

    /* Clamp D2: floor away from 0 (small-D2 is the quantisation regime
     * where 1/(d+D2) -> 1/d but D2/Q in the denominator blows the step
     * unless we keep D2 well above 0); cap at d. */
    double D2_eff = adapter_clamp(a->have_D2 ? a->D2_hat : (double)a->d,
                                   a->d2_floor, a->d2_ceil);

    /* (4) closed-form ell*. Until weights are seeded the step is a no-op. */
    double ell_raw;
    if (a->have_weights) {
        double B_cur = a->a_hat * (double)M_cur;
        double Q_cur = a->b_hat * (double)S_cur;
        ell_raw = adapter_optimal_ell(a->ell_cur, a->d, D2_eff, B_cur, Q_cur);
    } else {
        ell_raw = a->ell_cur;
    }

    /* (5) clamp ell to [ell_min, ell_max]. */
    double ell_star = adapter_clamp(ell_raw, a->ell_min, a->ell_max);

    /* (6) hysteresis. */
    int resized = adapter_should_resize(ell_star, a->ell_cur, a->eps);
    if (resized) a->ell_cur = ell_star;

    a->last_ell_raw    = ell_raw;
    a->last_ell_next   = ell_star;
    a->last_resized    = resized;
    a->last_overhead_s = now_s() - t_start;
    a->frame++;
    return a->ell_cur;
}

double adapter_blind_step(adapter_t *a,
                          int n, const double *px, const double *py, const double *pz) {
    double t_start = now_s();

    if ((a->frame % (int64_t)a->blind_period) != 0 || n <= 0) {
        a->last_ell_raw    = a->ell_cur;
        a->last_ell_next   = a->ell_cur;
        a->last_resized    = 0;
        a->last_overhead_s = now_s() - t_start;
        a->frame++;
        return a->ell_cur;
    }

    int K = a->blind_K > 0 ? a->blind_K : 9;
    if (K > 16) K = 16;
    double cands[16];
    int mid = K / 2;
    for (int k = 0; k < K; k++) {
        double mult = pow(a->blind_gamma, (double)(k - mid));
        double e = a->ell_cur * mult;
        if (e < a->ell_min) e = a->ell_min;
        if (e > a->ell_max) e = a->ell_max;
        cands[k] = e;
    }

    int best = mid;
    double best_t = 1e30;
    int pairs_cap = 16 * n;
    if (pairs_cap < 1024) pairs_cap = 1024;
    pair_t *pairs = (pair_t *)malloc(sizeof(pair_t) * (size_t)pairs_cap);
    if (!pairs) { a->last_overhead_s = now_s() - t_start; a->frame++; return a->ell_cur; }

    for (int k = 0; k < K; k++) {
        grid_t g;
        if (a->d == 3) grid_init3d(&g, n, a->domain_w, a->domain_h, a->domain_d, cands[k]);
        else           grid_init(&g, n, a->domain_w, a->domain_h, cands[k]);

        struct timespec ta, td;
        clock_gettime(CLOCK_MONOTONIC, &ta);
        if (a->d == 3) grid_build3d(&g, px, py, pz);
        else           grid_build(&g, px, py);
        int n_cand = 0;
        while (grid_broad_phase(&g, pairs, pairs_cap, &n_cand) != 0) {
            pairs_cap *= 2;
            pairs = (pair_t *)realloc(pairs, sizeof(pair_t) * (size_t)pairs_cap);
            if (!pairs) { grid_free(&g); a->last_overhead_s = now_s() - t_start; a->frame++; return a->ell_cur; }
        }
        int n_overlap;
        if (a->d == 3)
            narrow_phase3d(pairs, n_cand, px, py, pz, a->blind_radius, pairs, &n_overlap);
        else
            narrow_phase(pairs, n_cand, px, py, a->blind_radius, pairs, &n_overlap);
        clock_gettime(CLOCK_MONOTONIC, &td);

        double t_total = (double)(td.tv_sec - ta.tv_sec) +
                         (double)(td.tv_nsec - ta.tv_nsec) * 1e-9;
        if (t_total < best_t) { best_t = t_total; best = k; }
        grid_free(&g);
    }
    free(pairs);

    double ell_star = cands[best];
    int resized = adapter_should_resize(ell_star, a->ell_cur, a->eps);
    if (resized) a->ell_cur = ell_star;

    a->last_ell_raw    = ell_star;
    a->last_ell_next   = ell_star;
    a->last_resized    = resized;
    a->last_overhead_s = now_s() - t_start;
    a->frame++;
    return a->ell_cur;
}
