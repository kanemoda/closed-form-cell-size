/*
 * Phase 5 comparison harness (spec §9 -- the paper's results in miniature).
 *
 * For each scenario, at a MODERATE N (default 8K), run every method:
 *   - Ericson (ell = 2r)
 *   - Density / analytical (sqrt(V/N), clamped)
 *   - Static oracle (best single fixed ell, swept exhaustively)
 *   - Per-frame oracle (best ell PER frame, swept exhaustively)
 *   - Blind candidate search (the 9-candidate baseline)
 *   - Closed-form Mode A (D2 := d)
 *   - Closed-form Mode B (probed D2)
 *
 * For each method, report mean / p95 / p99 / max t_detect (us), plus the
 * mean-fraction gap to the per-frame oracle. Also report the closed-form
 * controller's measured per-frame overhead.
 *
 * Build:
 *   gcc -O2 -march=native -Wall -Wextra -Wno-unused-parameter -std=c11 -Isrc \
 *       -o compare_phase5 tools/compare_phase5.c \
 *       src/rng.c src/json.c src/config.c src/grid.c src/sim.c src/log.c \
 *       src/scenarios.c src/oracle.c src/baselines.c src/adapter.c -lm
 */

#define _POSIX_C_SOURCE 199309L

#include "config.h"
#include "sim.h"
#include "scenarios.h"
#include "grid.h"
#include "oracle.h"
#include "baselines.h"
#include "adapter.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    return da < db ? -1 : (da > db ? 1 : 0);
}
static double pct_sorted(const double *sorted, int n, double q) {
    if (n <= 0) return 0;
    int k = (int)floor(q * (double)(n - 1) + 0.5);
    if (k < 0) k = 0;
    if (k >= n) k = n - 1;
    return sorted[k];
}

typedef struct {
    const char *name;
    double mean_us, p95_us, p99_us, max_us;
    double overhead_us;    /* controller per-frame overhead (us); 0 for fixed/oracle */
    double mean_ell;       /* mean of chosen ell over the run */
} method_result_t;

/* Run sim with a FIXED cell size; return per-frame t_detect series in 'out',
 * length = num_frames. Allocated by caller. */
static void run_fixed(const config_t *cfg, double ell, double *t_detect_us,
                       double *mean_ell_out) {
    config_t mod = *cfg;
    mod.cell_size       = ell;
    mod.log_enabled     = 0;
    mod.log_csv_path[0] = '\0';
    mod.snapshot_csv_path[0] = '\0';
    mod.oracle_csv_path[0]   = '\0';
    mod.adapter_mode[0]      = '\0';

    sim_t s; sim_init(&s, &mod);
    scenario_init(&s, scenario_from_name(mod.scenario));
    for (int t = 0; t < mod.num_frames; t++) {
        sim_metrics_t m = { .frame = t };
        sim_step(&s, &m);
        t_detect_us[t] = (m.t_grid + m.t_broad + m.t_narrow) * 1e6;
    }
    sim_free(&s);
    if (mean_ell_out) *mean_ell_out = ell;
}

/* Run sim with the ADAPTIVE controller. Per-frame t_detect series in
 * t_detect_us. Controller's per-frame overhead in overhead_us (mean).
 * mean_ell_out = mean of ell_cur over the run.
 *
 * If 'fit' is non-NULL and well-formed (and the mode uses the closed form),
 * install its regression-fitted a, b as the controller's cost coefficients
 * instead of the per-frame sub-timer EMA. */
static void run_adapter(const config_t *cfg, adapter_mode_t mode,
                         const cost_fit_t *fit,
                         double *t_detect_us, double *overhead_us_mean,
                         double *mean_ell_out, double *ell_series) {
    config_t mod = *cfg;
    mod.log_enabled     = 0;
    mod.log_csv_path[0] = '\0';
    mod.snapshot_csv_path[0] = '\0';
    mod.oracle_csv_path[0]   = '\0';

    sim_t s; sim_init(&s, &mod);
    scenario_init(&s, scenario_from_name(mod.scenario));

    adapter_t a;
    adapter_init_default(&a, mode, 2, mod.domain_w, mod.domain_h,
                         mod.cell_size, mod.radius);

    /* Install regression-calibrated coefficients for the closed-form modes. */
    if (fit && fit->ok && fit->a > 0.0 && fit->b > 0.0 &&
        (mode == ADAPTER_MODE_A || mode == ADAPTER_MODE_B)) {
        adapter_set_calibrated_weights(&a, fit->a, fit->b);
    }

    double sum_oh = 0;
    double sum_ell = 0;
    for (int t = 0; t < mod.num_frames; t++) {
        sim_metrics_t m = { .frame = t };
        sim_step(&s, &m);
        t_detect_us[t] = (m.t_grid + m.t_broad + m.t_narrow) * 1e6;

        /* Spec §7 split: t_M includes broad-phase outer cell-iter;
         * t_S = inner pair-emit + narrow filter. */
        double t_M = m.t_M + m.t_broad_iter;
        double t_S = (m.t_broad - m.t_broad_iter) + m.t_narrow;
        if (t_S < 0) t_S = 0;
        double ell_next;
        if (mode == ADAPTER_MODE_BLIND)
            ell_next = adapter_blind_step(&a, s.n, s.px, s.py);
        else
            ell_next = adapter_step(&a, t_M, t_S, m.S, m.num_cells,
                                    s.n, s.px, s.py);
        sum_oh += a.last_overhead_s;
        sum_ell += a.ell_cur;
        if (ell_series) ell_series[t] = a.ell_cur;
        if (a.last_resized) grid_resize(&s.grid, ell_next);
    }
    sim_free(&s);
    if (overhead_us_mean) *overhead_us_mean = (sum_oh / mod.num_frames) * 1e6;
    if (mean_ell_out)     *mean_ell_out     = sum_ell / mod.num_frames;
}

/* Per-frame oracle: each frame we evaluate the full sweep and take argmin. */
static void run_per_frame_oracle(const config_t *cfg, const double *ells, int n_ells,
                                  double *t_detect_us, double *mean_ell_out) {
    config_t mod = *cfg;
    mod.log_enabled     = 0;
    mod.log_csv_path[0] = '\0';
    mod.snapshot_csv_path[0] = '\0';
    mod.oracle_csv_path[0]   = '\0';
    mod.adapter_mode[0]      = '\0';

    sim_t s; sim_init(&s, &mod);
    scenario_init(&s, scenario_from_name(mod.scenario));
    double sum_ell = 0;
    oracle_point_t *pts = (oracle_point_t *)malloc(sizeof(*pts) * n_ells);
    for (int t = 0; t < mod.num_frames; t++) {
        /* Note: per-frame oracle does NOT advance the sim at the best ell;
         * we measure t_detect from the oracle's own grid build (the spec
         * defines the per-frame oracle as the upper bound on any adaptive
         * method, evaluated non-perturbatively on each frame's snapshot).
         * So we still step the sim at the FIXED cfg.cell_size, then probe. */
        sim_step(&s, NULL);
        oracle_eval_sweep(s.n, s.px, s.py, mod.radius, mod.domain_w, mod.domain_h,
                          ells, n_ells, 3, pts);
        int best = 0;
        for (int i = 1; i < n_ells; i++)
            if (pts[i].t_detect < pts[best].t_detect) best = i;
        t_detect_us[t] = pts[best].t_detect * 1e6;
        sum_ell += pts[best].ell;
    }
    free(pts);
    sim_free(&s);
    if (mean_ell_out) *mean_ell_out = sum_ell / mod.num_frames;
}

/* One-time cost-model calibration (spec §2 cost model / §8 Pillar i).
 *
 * Step a fresh sim (same seed/config -- the trajectory is cell-size-
 * independent for ell >= 2r, since the overlap set the physics resolves does
 * not depend on the grid) to a few representative INTERIOR frames. At each,
 * sweep the candidate ells with the non-perturbing oracle and collect
 * (M, S, t_detect). Pool all samples and least-squares fit T ~= c0 + a*M +
 * b*S. The fitted a, b are the hardware constants the closed form needs; the
 * R^2 / residuals validate (or honestly refute) cost-model linearity.
 *
 * Note M = num_cells depends only on ell, so across frames M repeats at each
 * ell while S (clustering) varies -- that within-ell-across-frame spread in S
 * is what makes the two-predictor fit identifiable. */
static cost_fit_t calibrate_cost_model(const config_t *cfg,
                                       const double *ells, int n_ells,
                                       int n_calib_frames, int K) {
    config_t mod = *cfg;
    mod.log_enabled          = 0;
    mod.log_csv_path[0]      = '\0';
    mod.snapshot_csv_path[0] = '\0';
    mod.oracle_csv_path[0]   = '\0';
    mod.adapter_mode[0]      = '\0';

    sim_t s; sim_init(&s, &mod);
    scenario_init(&s, scenario_from_name(mod.scenario));

    int F = mod.num_frames;
    if (n_calib_frames < 1)  n_calib_frames = 1;
    if (n_calib_frames > 64) n_calib_frames = 64;

    /* Evenly spaced interior target frames: skip frame 0's placement
     * transient and the last frame. */
    int targets[64];
    for (int k = 0; k < n_calib_frames; k++) {
        int tf = (int)((double)(k + 1) * (double)F / (double)(n_calib_frames + 1) + 0.5);
        if (tf < 0)      tf = 0;
        if (tf > F - 1)  tf = F - 1;
        targets[k] = tf;
    }

    int cap = n_calib_frames * n_ells;
    double *Ms = (double *)malloc(sizeof(double) * (size_t)cap);
    double *Ss = (double *)malloc(sizeof(double) * (size_t)cap);
    double *Ts = (double *)malloc(sizeof(double) * (size_t)cap);
    oracle_point_t *pts = (oracle_point_t *)malloc(sizeof(*pts) * (size_t)n_ells);
    int ns = 0, tk = 0;

    for (int t = 0; t < F && tk < n_calib_frames; t++) {
        sim_step(&s, NULL);   /* advance physics at the fixed cfg cell size */
        if (t != targets[tk]) continue;
        oracle_eval_sweep(s.n, s.px, s.py, mod.radius, mod.domain_w, mod.domain_h,
                          ells, n_ells, K, pts);
        for (int i = 0; i < n_ells && ns < cap; i++) {
            Ms[ns] = (double)pts[i].num_cells;
            Ss[ns] = (double)pts[i].S;
            Ts[ns] = pts[i].t_detect;
            ns++;
        }
        tk++;
    }

    cost_fit_t fit = adapter_fit_cost_model(Ms, Ss, Ts, ns);
    free(Ms); free(Ss); free(Ts); free(pts);
    sim_free(&s);
    return fit;
}

/* Static oracle: try each candidate ell as a full fixed run; pick the one
 * with lowest mean t_detect. Then return its per-frame series. */
static double find_static_best_ell(const config_t *cfg, const double *ells, int n_ells) {
    double best_mean = 1e30;
    double best_ell  = ells[0];
    double *buf = (double *)malloc(sizeof(double) * cfg->num_frames);
    for (int i = 0; i < n_ells; i++) {
        double ell_mean;
        run_fixed(cfg, ells[i], buf, &ell_mean);
        double sum = 0;
        for (int t = 0; t < cfg->num_frames; t++) sum += buf[t];
        double mean = sum / cfg->num_frames;
        if (mean < best_mean) { best_mean = mean; best_ell = ells[i]; }
    }
    free(buf);
    return best_ell;
}

static void stats(const double *series, int n, method_result_t *r) {
    double *s = (double *)malloc(sizeof(double) * n);
    memcpy(s, series, sizeof(double) * n);
    double sum = 0;
    for (int i = 0; i < n; i++) sum += s[i];
    r->mean_us = sum / n;
    qsort(s, n, sizeof(double), cmp_double);
    r->p95_us = pct_sorted(s, n, 0.95);
    r->p99_us = pct_sorted(s, n, 0.99);
    r->max_us = s[n - 1];
    free(s);
}

static void run_scenario(scenario_id_t scen_id, int N, int frames) {
    config_t cfg;
    config_default(&cfg);
    cfg.N           = N;
    cfg.domain_w    = 800.0;
    cfg.domain_h    = 600.0;
    cfg.radius      = 0.5;
    cfg.dt          = 1.0/60.0;
    cfg.init_speed  = 10.0;
    cfg.restitution = 0.95;
    cfg.num_frames  = frames;
    cfg.seed        = 31;
    cfg.cell_size   = baseline_uniform_analytical(N, cfg.domain_w, cfg.domain_h, cfg.radius);
    snprintf(cfg.scenario, sizeof(cfg.scenario), "%s", scenario_name(scen_id));

    fprintf(stderr,
        "\n=== %s, N=%d, frames=%d ===\n",
        scenario_name(scen_id), N, frames);
    fflush(stderr);

    /* Candidate ell sweep for static & per-frame oracles. */
    double ells[32];
    int n_ells = oracle_make_geo_sweep(2.0 * cfg.radius, 20.0, 1.25, ells, 32);

    /* One-time regression calibration of the closed form's cost coefficients
     * (spec §2/§8). Fit over the same candidate sweep across a few
     * representative frames; the fitted a, b feed Closed-A and Closed-B. */
    const int n_calib_frames = 8, calib_K = 5;
    fprintf(stderr, "  Calibrating cost model (%d frames x %d ells, K=%d) ...\n",
            n_calib_frames, n_ells, calib_K);
    fflush(stderr);
    cost_fit_t fit = calibrate_cost_model(&cfg, ells, n_ells, n_calib_frames, calib_K);
    if (fit.ok) {
        fprintf(stderr,
            "  cost model  T = c0 + a*M + b*S :  a=%.3e s/cell  b=%.3e s/unitS  "
            "c0=%.3e s   [n=%d]\n",
            fit.a, fit.b, fit.c0, fit.n);
        fprintf(stderr,
            "    Pillar-i fit: R^2=%.4f  RMSE=%.2f us  mean|resid|=%.2f us (%.1f%%)  "
            "max|resid|=%.2f us\n",
            fit.r2, fit.rmse * 1e6, fit.mean_abs_resid * 1e6,
            fit.mean_rel_resid * 100.0, fit.max_abs_resid * 1e6);
        if (fit.a <= 0.0 || fit.b <= 0.0)
            fprintf(stderr,
                "    WARNING: non-positive coefficient -- closed form falls back "
                "to per-frame sub-timer EMA for this scenario.\n");
    } else {
        fprintf(stderr, "  cost-model calibration FAILED (degenerate fit); "
                        "closed form falls back to per-frame sub-timer EMA.\n");
    }
    const cost_fit_t *fitp = (fit.ok && fit.a > 0.0 && fit.b > 0.0) ? &fit : NULL;
    fflush(stderr);

    /* Per-method t_detect series (us) */
    double *t = (double *)malloc(sizeof(double) * frames);
    method_result_t res[8] = {0};
    int n_methods = 0;

    fprintf(stderr, "  Ericson ...\n");      fflush(stderr);
    double ell_eric = baseline_ericson(cfg.radius);
    run_fixed(&cfg, ell_eric, t, &res[n_methods].mean_ell);
    stats(t, frames, &res[n_methods]); res[n_methods].name = "Ericson(2r)";
    n_methods++;

    fprintf(stderr, "  Density baseline ...\n"); fflush(stderr);
    double ell_dens = baseline_uniform_analytical(cfg.N, cfg.domain_w, cfg.domain_h, cfg.radius);
    run_fixed(&cfg, ell_dens, t, &res[n_methods].mean_ell);
    stats(t, frames, &res[n_methods]); res[n_methods].name = "Density(sqrt(V/N))";
    n_methods++;

    fprintf(stderr, "  Static oracle (sweep) ...\n"); fflush(stderr);
    double ell_stat = find_static_best_ell(&cfg, ells, n_ells);
    run_fixed(&cfg, ell_stat, t, &res[n_methods].mean_ell);
    stats(t, frames, &res[n_methods]); res[n_methods].name = "Static-oracle";
    n_methods++;

    fprintf(stderr, "  Blind 9-candidate ...\n"); fflush(stderr);
    run_adapter(&cfg, ADAPTER_MODE_BLIND, NULL, t, &res[n_methods].overhead_us,
                &res[n_methods].mean_ell, NULL);
    stats(t, frames, &res[n_methods]); res[n_methods].name = "Blind-9";
    n_methods++;

    fprintf(stderr, "  Closed-form Mode A (calibrated) ...\n"); fflush(stderr);
    run_adapter(&cfg, ADAPTER_MODE_A, fitp, t, &res[n_methods].overhead_us,
                &res[n_methods].mean_ell, NULL);
    stats(t, frames, &res[n_methods]); res[n_methods].name = "Closed-A";
    n_methods++;

    fprintf(stderr, "  Closed-form Mode B (calibrated) ...\n"); fflush(stderr);
    run_adapter(&cfg, ADAPTER_MODE_B, fitp, t, &res[n_methods].overhead_us,
                &res[n_methods].mean_ell, NULL);
    stats(t, frames, &res[n_methods]); res[n_methods].name = "Closed-B";
    n_methods++;

    fprintf(stderr, "  Per-frame oracle ...\n"); fflush(stderr);
    run_per_frame_oracle(&cfg, ells, n_ells, t, &res[n_methods].mean_ell);
    stats(t, frames, &res[n_methods]); res[n_methods].name = "PF-oracle";
    n_methods++;

    free(t);

    /* Look up PF-oracle stats for ratio columns. */
    method_result_t pf = res[n_methods - 1];

    /* Print results table. */
    fprintf(stderr,
        "\n  method            |  mean_us | p95_us | p99_us | max_us |"
        "  mean_ell | oh_us | mean/PF\n");
    fprintf(stderr,
        "  ------------------+----------+--------+--------+--------+"
        "----------+-------+--------\n");
    for (int i = 0; i < n_methods; i++) {
        fprintf(stderr,
            "  %-17s |  %7.2f | %6.2f | %6.2f | %6.2f |  %7.3f | %5.2f |  %.3f\n",
            res[i].name,
            res[i].mean_us, res[i].p95_us, res[i].p99_us, res[i].max_us,
            res[i].mean_ell, res[i].overhead_us,
            res[i].mean_us / pf.mean_us);
    }
}

int main(int argc, char **argv) {
    int N = 8000;
    int frames = 200;
    if (argc > 1) N = atoi(argv[1]);
    if (argc > 2) frames = atoi(argv[2]);

    scenario_id_t scenarios[] = {
        SCEN_STABLE, SCEN_COLLAPSE, SCEN_EXPLOSION, SCEN_VORTEX,
        SCEN_FUNNEL, SCEN_PULSE, SCEN_MULTICLUSTER, SCEN_STREAM
    };
    int n_scen = (int)(sizeof(scenarios)/sizeof(scenarios[0]));

    fprintf(stderr, "===================================================================\n");
    fprintf(stderr, "Phase 5 method comparison (N=%d, frames=%d, 800x600 domain)\n", N, frames);
    fprintf(stderr, "===================================================================\n");

    /* multicluster's initial placement may exhaust at high N in narrow disks --
     * fall back to a more modest N for that scenario, mirroring diag.c's note. */
    for (int i = 0; i < n_scen; i++) {
        int N_use = N;
        if (scenarios[i] == SCEN_MULTICLUSTER && N > 6000) N_use = 6000;
        run_scenario(scenarios[i], N_use, frames);
    }
    return 0;
}
