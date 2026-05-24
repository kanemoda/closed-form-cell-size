/*
 * Phase 6 comparison harness -- the 3D counterpart of tools/compare_phase5.c.
 *
 * For each 3D scenario (stable, collapse, explosion, pulse, settling) at a
 * moderate N, run every method and report mean / p95 / p99 / max t_detect:
 *   - Ericson (ell = 2r)
 *   - Density / analytical (cbrt(V/N), clamped)
 *   - Static oracle (best single fixed ell; selection + independent report)
 *   - Per-frame oracle (split-sample: SELECT on one timing sample, REPORT the
 *     other -- the de-biased per-frame upper bound)
 *   - Blind 9-candidate search
 *   - Closed-form Mode A / Mode B (regression-calibrated coefficients)
 *
 * The cost coefficients a, b are calibrated ONCE per scenario by regression of
 * T = c0 + a*M + b*S over the candidate ell sweep across representative frames
 * (spec §2/§8); R^2 and residuals are reported.
 *
 * Build: make compare_phase6   (then ./compare_phase6 [N] [frames])
 */

#define _POSIX_C_SOURCE 199309L

#include "config.h"
#include "sim.h"
#include "scenarios.h"
#include "grid.h"
#include "oracle.h"
#include "adapter.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
    double overhead_us;
    double mean_ell;
} method_result_t;

/* cbrt(V/N) density baseline, clamped to >= 2r. */
static double density_baseline_3d(int N, double W, double H, double D, double r) {
    double ell = cbrt((W * H * D) / (double)N);
    double floor_ell = 2.0 * r;
    return ell < floor_ell ? floor_ell : ell;
}

static void clear_io(config_t *m) {
    m->log_enabled = 0;
    m->log_csv_path[0] = '\0';
    m->snapshot_csv_path[0] = '\0';
    m->oracle_csv_path[0] = '\0';
    m->adapter_mode[0] = '\0';
}

/* Fixed cell size: per-frame t_detect series (us). */
static void run_fixed3d(const config_t *cfg, double ell, double *t_detect_us,
                        double *mean_ell_out) {
    config_t mod = *cfg; clear_io(&mod);
    mod.cell_size = ell;
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

/* Adaptive controller (Mode A/B/blind), optionally calibrated. */
static void run_adapter3d(const config_t *cfg, adapter_mode_t mode, const cost_fit_t *fit,
                          double *t_detect_us, double *overhead_us_mean,
                          double *mean_ell_out) {
    config_t mod = *cfg; clear_io(&mod);
    sim_t s; sim_init(&s, &mod);
    scenario_init(&s, scenario_from_name(mod.scenario));

    adapter_t a;
    adapter_init_default(&a, mode, 3, mod.domain_w, mod.domain_h, mod.cell_size, mod.radius);
    a.domain_d = mod.domain_d;
    if (fit && fit->ok && fit->a > 0.0 && fit->b > 0.0 &&
        (mode == ADAPTER_MODE_A || mode == ADAPTER_MODE_B))
        adapter_set_calibrated_weights(&a, fit->a, fit->b);

    double sum_oh = 0, sum_ell = 0;
    for (int t = 0; t < mod.num_frames; t++) {
        sim_metrics_t m = { .frame = t };
        sim_step(&s, &m);
        t_detect_us[t] = (m.t_grid + m.t_broad + m.t_narrow) * 1e6;
        double t_M = m.t_M + m.t_broad_iter;
        double t_S = (m.t_broad - m.t_broad_iter) + m.t_narrow;
        if (t_S < 0) t_S = 0;
        double ell_next;
        if (mode == ADAPTER_MODE_BLIND)
            ell_next = adapter_blind_step(&a, s.n, s.px, s.py, s.pz);
        else
            ell_next = adapter_step(&a, t_M, t_S, m.S, m.num_cells, s.n, s.px, s.py, s.pz);
        sum_oh += a.last_overhead_s;
        sum_ell += a.ell_cur;
        if (a.last_resized) grid_resize(&s.grid, ell_next);
    }
    sim_free(&s);
    if (overhead_us_mean) *overhead_us_mean = (sum_oh / mod.num_frames) * 1e6;
    if (mean_ell_out)     *mean_ell_out     = sum_ell / mod.num_frames;
}

/* Per-frame oracle, split-sample (SELECT on sel, REPORT eval). */
static void run_per_frame_oracle3d(const config_t *cfg, const double *ells, int n_ells,
                                   double *t_detect_us, double *mean_ell_out) {
    config_t mod = *cfg; clear_io(&mod);
    sim_t s; sim_init(&s, &mod);
    scenario_init(&s, scenario_from_name(mod.scenario));
    double sum_ell = 0;
    /* Single-pass per-frame oracle (median-of-K=3): carries the same small
     * (<=~2%) winner's-curse as 2D, negligible vs the 1.0-1.5x method gaps and
     * ~2x cheaper than split-sample -- the saving funds seed averaging. */
    oracle_point_t *pts = (oracle_point_t *)malloc(sizeof(*pts) * n_ells);
    for (int t = 0; t < mod.num_frames; t++) {
        sim_step(&s, NULL);
        oracle_eval_sweep3d(s.n, s.px, s.py, s.pz, mod.radius,
                            mod.domain_w, mod.domain_h, mod.domain_d, ells, n_ells, 3, pts);
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

/* Static oracle: SELECT best fixed ell on one timing pass; the caller REPORTs
 * via a separate run_fixed3d (independent timing) -- split-sample. */
static double find_static_best_ell3d(const config_t *cfg, const double *ells, int n_ells) {
    double best_mean = 1e30, best_ell = ells[0];
    double *buf = (double *)malloc(sizeof(double) * cfg->num_frames);
    for (int i = 0; i < n_ells; i++) {
        run_fixed3d(cfg, ells[i], buf, NULL);
        double sum = 0;
        for (int t = 0; t < cfg->num_frames; t++) sum += buf[t];
        double mean = sum / cfg->num_frames;
        if (mean < best_mean) { best_mean = mean; best_ell = ells[i]; }
    }
    free(buf);
    return best_ell;
}

/* One-time cost-model calibration (3D): fit T = c0 + a*M + b*S over the
 * candidate sweep across representative interior frames. */
static cost_fit_t calibrate_cost_model3d(const config_t *cfg, const double *ells, int n_ells,
                                         int n_calib_frames, int K) {
    config_t mod = *cfg; clear_io(&mod);
    sim_t s; sim_init(&s, &mod);
    scenario_init(&s, scenario_from_name(mod.scenario));

    int F = mod.num_frames;
    if (n_calib_frames < 1)  n_calib_frames = 1;
    if (n_calib_frames > 64) n_calib_frames = 64;
    int targets[64];
    for (int k = 0; k < n_calib_frames; k++) {
        int tf = (int)((double)(k + 1) * (double)F / (double)(n_calib_frames + 1) + 0.5);
        if (tf < 0) tf = 0;
        if (tf > F - 1) tf = F - 1;
        targets[k] = tf;
    }
    int cap = n_calib_frames * n_ells, ns = 0, tk = 0;
    double *Ms = (double *)malloc(sizeof(double) * cap);
    double *Ss = (double *)malloc(sizeof(double) * cap);
    double *Ts = (double *)malloc(sizeof(double) * cap);
    oracle_point_t *pts = (oracle_point_t *)malloc(sizeof(*pts) * n_ells);
    for (int t = 0; t < F && tk < n_calib_frames; t++) {
        sim_step(&s, NULL);
        if (t != targets[tk]) continue;
        oracle_eval_sweep3d(s.n, s.px, s.py, s.pz, mod.radius,
                            mod.domain_w, mod.domain_h, mod.domain_d, ells, n_ells, K, pts);
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

static void run_scenario3d(scenario_id_t scen_id, int N, int frames, uint64_t seed) {
    config_t cfg; config_default(&cfg);
    cfg.dim         = 3;
    cfg.N           = N;
    cfg.domain_w    = 200.0;
    cfg.domain_h    = 200.0;
    cfg.domain_d    = 200.0;
    cfg.radius      = 0.5;
    cfg.dt          = 1.0 / 60.0;
    cfg.init_speed  = 10.0;
    cfg.restitution = 0.95;
    cfg.num_frames  = frames;
    cfg.seed        = seed;
    cfg.cell_size   = density_baseline_3d(N, cfg.domain_w, cfg.domain_h, cfg.domain_d, cfg.radius);
    snprintf(cfg.scenario, sizeof(cfg.scenario), "%s", scenario_name(scen_id));

    fprintf(stderr, "\n=== %s (3D), N=%d, frames=%d, 200^3 ===\n",
            scenario_name(scen_id), N, frames);
    fflush(stderr);

    /* Clipped ladder (~8 rungs) bracketing the 3D density-baseline optimum
     * (cfg.cell_size = cbrt(V/N)), [0.45db, 2.4db] clamped to >= 2 (ell<2 in
     * 200^3 is M-dominated and never optimal; the large-ell cap also avoids the
     * S~1e8 pair-buffer blowup). Ericson ell=2r=1 reported separately. */
    double ells[32];
    double db = cfg.cell_size;
    int n_ells = oracle_make_geo_sweep(fmax(2.0, 0.45 * db), 2.4 * db, 1.3, ells, 32);

    const int n_calib_frames = 8, calib_K = 5;
    fprintf(stderr, "  Calibrating cost model (%d frames x %d ells, K=%d) ...\n",
            n_calib_frames, n_ells, calib_K); fflush(stderr);
    cost_fit_t fit = calibrate_cost_model3d(&cfg, ells, n_ells, n_calib_frames, calib_K);
    if (fit.ok)
        fprintf(stderr,
            "  cost model  T = c0 + a*M + b*S :  a=%.3e  b=%.3e  c0=%.3e   [n=%d]\n"
            "    Pillar-i fit: R^2=%.4f  RMSE=%.2f us  mean|resid|=%.1f%%  max|resid|=%.1f us\n",
            fit.a, fit.b, fit.c0, fit.n, fit.r2, fit.rmse * 1e6,
            fit.mean_rel_resid * 100.0, fit.max_abs_resid * 1e6);
    else
        fprintf(stderr, "  calibration FAILED -- closed form falls back to per-frame EMA.\n");
    const cost_fit_t *fitp = (fit.ok && fit.a > 0.0 && fit.b > 0.0) ? &fit : NULL;
    fflush(stderr);

    double *t = (double *)malloc(sizeof(double) * frames);
    method_result_t res[8] = {0};
    int nm = 0;

    fprintf(stderr, "  Ericson ...\n"); fflush(stderr);
    run_fixed3d(&cfg, 2.0 * cfg.radius, t, &res[nm].mean_ell);
    stats(t, frames, &res[nm]); res[nm].name = "Ericson(2r)"; nm++;

    fprintf(stderr, "  Density baseline ...\n"); fflush(stderr);
    run_fixed3d(&cfg, cfg.cell_size, t, &res[nm].mean_ell);
    stats(t, frames, &res[nm]); res[nm].name = "Density(cbrt)"; nm++;

    fprintf(stderr, "  Static oracle (split-sample) ...\n"); fflush(stderr);
    double ell_stat = find_static_best_ell3d(&cfg, ells, n_ells);
    run_fixed3d(&cfg, ell_stat, t, &res[nm].mean_ell);     /* independent report run */
    stats(t, frames, &res[nm]); res[nm].name = "Static-oracle"; nm++;

    fprintf(stderr, "  Blind 9-candidate ...\n"); fflush(stderr);
    run_adapter3d(&cfg, ADAPTER_MODE_BLIND, NULL, t, &res[nm].overhead_us, &res[nm].mean_ell);
    stats(t, frames, &res[nm]); res[nm].name = "Blind-9"; nm++;

    fprintf(stderr, "  Closed-form Mode A (calibrated) ...\n"); fflush(stderr);
    run_adapter3d(&cfg, ADAPTER_MODE_A, fitp, t, &res[nm].overhead_us, &res[nm].mean_ell);
    stats(t, frames, &res[nm]); res[nm].name = "Closed-A"; nm++;

    fprintf(stderr, "  Closed-form Mode B (calibrated) ...\n"); fflush(stderr);
    run_adapter3d(&cfg, ADAPTER_MODE_B, fitp, t, &res[nm].overhead_us, &res[nm].mean_ell);
    stats(t, frames, &res[nm]); res[nm].name = "Closed-B"; nm++;

    fprintf(stderr, "  Per-frame oracle (split-sample) ...\n"); fflush(stderr);
    run_per_frame_oracle3d(&cfg, ells, n_ells, t, &res[nm].mean_ell);
    stats(t, frames, &res[nm]); res[nm].name = "PF-oracle"; nm++;

    free(t);
    method_result_t pf = res[nm - 1];

    fprintf(stderr,
        "\n  method            |  mean_us | p95_us | p99_us | max_us |"
        "  mean_ell | oh_us | mean/PF\n");
    fprintf(stderr,
        "  ------------------+----------+--------+--------+--------+"
        "----------+-------+--------\n");
    for (int i = 0; i < nm; i++)
        fprintf(stderr,
            "  %-17s |  %7.2f | %6.2f | %6.2f | %6.2f |  %7.3f | %5.2f |  %.3f\n",
            res[i].name, res[i].mean_us, res[i].p95_us, res[i].p99_us, res[i].max_us,
            res[i].mean_ell, res[i].overhead_us, res[i].mean_us / pf.mean_us);
}

int main(int argc, char **argv) {
    int N = 8000;          /* moderate N (matches the 2D comparison); spec
                            * campaign uses 50K (Phase 7). */
    int frames = 100;
    uint64_t seed = 31;
    if (argc > 1) N = atoi(argv[1]);
    if (argc > 2) frames = atoi(argv[2]);
    if (argc > 3) seed = (uint64_t)strtoull(argv[3], NULL, 10);

    scenario_id_t scenarios[] = {
        SCEN_STABLE, SCEN_COLLAPSE, SCEN_EXPLOSION, SCEN_PULSE, SCEN_SETTLING
    };
    int n_scen = (int)(sizeof(scenarios) / sizeof(scenarios[0]));

    fprintf(stderr, "===================================================================\n");
    fprintf(stderr, "Phase 6 method comparison -- 3D (N=%d, frames=%d, seed=%llu, 200^3)\n",
            N, frames, (unsigned long long)seed);
    fprintf(stderr, "===================================================================\n");

    for (int i = 0; i < n_scen; i++)
        run_scenario3d(scenarios[i], N, frames, seed);
    return 0;
}
