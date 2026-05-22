/*
 * Phase 4 verification gates (spec section 10):
 *   G1. Oracle determinism: same frozen snapshot at the same ell twice
 *       -> bit-identical S, num_cells, candidate_pairs.
 *   G2. Non-perturbing: a run with the per-frame oracle ACTIVE produces
 *       the same per-frame deterministic state (S, num_cells, candidate
 *       pairs, collisions, KE, momentum) as the same run without it.
 *   G3. The sweep brackets the optimum: for each scenario, on at least
 *       one representative frame the argmin of t_detect(ell) is interior
 *       to the sweep (not at ell_min or ell_max). If a scenario's
 *       optimum sits on a boundary, widen the sweep and re-run.
 *   G4. Static oracle beats Ericson (ell = 2r) on mean t_detect, for
 *       each scenario.
 *
 * Plus: prints per-scenario (ell, S) sweep tables on a few representative
 * frames as raw material for the S ~ ell^D2 power-law check (Phase 8).
 */

#include "config.h"
#include "sim.h"
#include "scenarios.h"
#include "grid.h"
#include "oracle.h"
#include "baselines.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ===================================================== shared test config */

static config_t std_cfg(scenario_id_t id) {
    config_t cfg;
    config_default(&cfg);
    cfg.N           = 300;
    cfg.domain_w    = 60.0;
    cfg.domain_h    = 60.0;
    cfg.radius      = 0.5;
    cfg.cell_size   = 2.0;
    cfg.dt          = 1.0 / 60.0;
    cfg.init_speed  = 10.0;
    cfg.restitution = 0.95;
    cfg.num_frames  = 60;
    cfg.seed        = 31;
    strncpy(cfg.scenario, scenario_name(id), sizeof(cfg.scenario) - 1);
    cfg.scenario[sizeof(cfg.scenario) - 1] = '\0';
    return cfg;
}

/* The default sweep used by gates G1, G3 and the (ell, S) tables.
 * Wide enough to bracket the optimum on every scenario at our test scale. */
static int make_sweep(const config_t *cfg, double *ells, int max_pts) {
    double ell_min = 2.0 * cfg->radius;
    double ell_max = 0.5 * fmin(cfg->domain_w, cfg->domain_h);
    return oracle_make_geo_sweep(ell_min, ell_max, 1.25, ells, max_pts);
}

/* =========================================================== G1: determinism */

static int g1_oracle_deterministic(void) {
    /* Build a non-trivial snapshot: run sim 20 frames, then freeze. */
    config_t cfg = std_cfg(SCEN_STABLE);
    cfg.num_frames = 20;
    sim_t s; sim_init(&s, &cfg); scenario_init(&s, SCEN_STABLE);
    for (int t = 0; t < cfg.num_frames; t++) sim_step(&s, NULL);

    int ok = 1;
    /* Three test ells: floor (2r), a typical middle value, and a wide one. */
    double test_ells[3] = { 2.0 * cfg.radius, 4.0, 8.0 };
    for (int k = 0; k < 3; k++) {
        oracle_point_t p1, p2;
        oracle_eval(s.n, s.px, s.py, cfg.radius,
                    cfg.domain_w, cfg.domain_h, test_ells[k], 5, &p1);
        oracle_eval(s.n, s.px, s.py, cfg.radius,
                    cfg.domain_w, cfg.domain_h, test_ells[k], 5, &p2);
        int eq = (p1.S == p2.S) && (p1.num_cells == p2.num_cells)
              && (p1.candidate_pairs == p2.candidate_pairs);
        if (!eq) ok = 0;
        fprintf(stderr,
            "  ell=%.3f : S=%ld/%ld  cells=%d/%d  cands=%d/%d  %s\n",
            test_ells[k], p1.S, p2.S, p1.num_cells, p2.num_cells,
            p1.candidate_pairs, p2.candidate_pairs, eq ? "OK" : "FAIL");
    }
    sim_free(&s);
    fprintf(stderr, "  G1 oracle determinism: %s\n", ok ? "OK" : "FAIL");
    return ok;
}

/* =========================================================== G2: non-perturbing */

static int g2_non_perturbing(void) {
    /* Run the same scenario twice with the same (seed, cfg). The "active"
     * variant runs oracle_eval_sweep after every sim_step but discards
     * the results. Compare deterministic per-frame metrics. */
    config_t cfg = std_cfg(SCEN_STABLE);
    cfg.num_frames = 30;

    double ells[32];
    int n_ells = make_sweep(&cfg, ells, 32);

    int frames = cfg.num_frames;
    long   *S_a   = malloc(sizeof(long)   * (size_t)frames);
    long   *S_b   = malloc(sizeof(long)   * (size_t)frames);
    int    *NC_a  = malloc(sizeof(int)    * (size_t)frames);
    int    *NC_b  = malloc(sizeof(int)    * (size_t)frames);
    int    *CP_a  = malloc(sizeof(int)    * (size_t)frames);
    int    *CP_b  = malloc(sizeof(int)    * (size_t)frames);
    int    *CL_a  = malloc(sizeof(int)    * (size_t)frames);
    int    *CL_b  = malloc(sizeof(int)    * (size_t)frames);
    double *KE_a  = malloc(sizeof(double) * (size_t)frames);
    double *KE_b  = malloc(sizeof(double) * (size_t)frames);
    double *PM_a  = malloc(sizeof(double) * (size_t)frames);
    double *PM_b  = malloc(sizeof(double) * (size_t)frames);

    /* Run A: no oracle. */
    {
        sim_t s; sim_init(&s, &cfg); scenario_init(&s, SCEN_STABLE);
        for (int t = 0; t < frames; t++) {
            sim_metrics_t m = { .frame = t };
            sim_step(&s, &m);
            S_a [t] = m.S;   NC_a[t] = m.num_cells;
            CP_a[t] = m.candidate_pairs; CL_a[t] = m.collisions;
            KE_a[t] = m.kinetic_energy;  PM_a[t] = m.momentum_magnitude;
        }
        sim_free(&s);
    }

    /* Run B: oracle eval after every step, results discarded. */
    {
        sim_t s; sim_init(&s, &cfg); scenario_init(&s, SCEN_STABLE);
        oracle_point_t pts[32];
        for (int t = 0; t < frames; t++) {
            sim_metrics_t m = { .frame = t };
            sim_step(&s, &m);
            S_b [t] = m.S;   NC_b[t] = m.num_cells;
            CP_b[t] = m.candidate_pairs; CL_b[t] = m.collisions;
            KE_b[t] = m.kinetic_energy;  PM_b[t] = m.momentum_magnitude;
            oracle_eval_sweep(s.n, s.px, s.py, cfg.radius,
                              cfg.domain_w, cfg.domain_h,
                              ells, n_ells, 2, pts);
        }
        sim_free(&s);
    }

    int mism = 0;
    for (int t = 0; t < frames; t++) {
        if (S_a[t]  != S_b[t])  mism++;
        if (NC_a[t] != NC_b[t]) mism++;
        if (CP_a[t] != CP_b[t]) mism++;
        if (CL_a[t] != CL_b[t]) mism++;
        if (KE_a[t] != KE_b[t]) mism++;   /* exact double match required */
        if (PM_a[t] != PM_b[t]) mism++;
    }
    int ok = (mism == 0);
    fprintf(stderr,
        "  G2 non-perturbing: %d deterministic mismatches over %d frames * 6 cols = %d  %s\n",
        mism, frames, frames * 6, ok ? "OK" : "FAIL");

    free(S_a); free(S_b); free(NC_a); free(NC_b); free(CP_a); free(CP_b);
    free(CL_a); free(CL_b); free(KE_a); free(KE_b); free(PM_a); free(PM_b);
    return ok;
}

/* =========================================================== G3: sweep brackets */

static int g3_sweep_brackets_one(scenario_id_t id) {
    config_t cfg = std_cfg(id);
    cfg.num_frames = 30;

    double ells[32];
    int n_ells = make_sweep(&cfg, ells, 32);

    sim_t s; sim_init(&s, &cfg); scenario_init(&s, id);
    oracle_point_t pts[32];
    int interior = 0, boundary = 0;
    int last_best = -1;
    for (int t = 0; t < cfg.num_frames; t++) {
        sim_step(&s, NULL);
        oracle_eval_sweep(s.n, s.px, s.py, cfg.radius,
                          cfg.domain_w, cfg.domain_h,
                          ells, n_ells, 5, pts);
        int best = 0;
        for (int i = 1; i < n_ells; i++) {
            if (pts[i].t_detect < pts[best].t_detect) best = i;
        }
        last_best = best;
        if (best == 0 || best == n_ells - 1) boundary++;
        else                                  interior++;
    }
    int ok = (interior > 0);   /* at least one frame's optimum is interior */
    fprintf(stderr,
        "  %-12s  n_ells=%d  interior_argmins=%d  boundary=%d  (last best idx=%d, ell=%.3f)  %s\n",
        scenario_name(id), n_ells, interior, boundary, last_best,
        last_best >= 0 ? ells[last_best] : -1.0, ok ? "OK" : "FAIL");
    sim_free(&s);
    return ok;
}

static int g3_sweep_brackets_all(void) {
    int ok = 1;
    for (int id = 0; id < SCEN_COUNT; id++) {
        ok &= g3_sweep_brackets_one((scenario_id_t)id);
    }
    return ok;
}

/* =========================================================== G4: oracle beats Ericson */

static int g4_oracle_beats_ericson_one(scenario_id_t id) {
    config_t cfg = std_cfg(id);
    cfg.num_frames = 40;

    double ells[32];
    int n_ells = make_sweep(&cfg, ells, 32);

    static_oracle_result_t res;
    static_oracle_run(&cfg, ells, n_ells, &res);

    sim_run_summary_t er;
    run_sim_with_cell_size(&cfg, baseline_ericson(cfg.radius), &er);

    double analytical = baseline_uniform_analytical(cfg.N, cfg.domain_w,
                                                    cfg.domain_h, cfg.radius);
    sim_run_summary_t an;
    run_sim_with_cell_size(&cfg, analytical, &an);

    int ok = (res.best_mean_detect <= er.mean_t_detect);
    fprintf(stderr,
        "  %-12s  oracle ell=%.3f  t=%.3f us  |  ericson(2r) t=%.3f us  |  "
        "analytical(%.3f) t=%.3f us  |  speedup=%.2fx  %s\n",
        scenario_name(id),
        res.best_ell,        res.best_mean_detect  * 1e6,
        er.mean_t_detect    * 1e6,
        analytical,         an.mean_t_detect       * 1e6,
        er.mean_t_detect / res.best_mean_detect,
        ok ? "OK" : "FAIL");

    static_oracle_free(&res);
    return ok;
}

static int g4_oracle_beats_ericson_all(void) {
    int ok = 1;
    for (int id = 0; id < SCEN_COUNT; id++) {
        ok &= g4_oracle_beats_ericson_one((scenario_id_t)id);
    }
    return ok;
}

/* =============================== (ell, S) sweep tables for the summary */

static void print_S_tables_one(scenario_id_t id) {
    config_t cfg = std_cfg(id);
    cfg.num_frames = 60;

    double ells[32];
    int n_ells = make_sweep(&cfg, ells, 32);

    /* Three representative frames -- early, mid, late. */
    int sample_frames[3] = { 5, cfg.num_frames / 2, cfg.num_frames - 1 };

    sim_t s; sim_init(&s, &cfg); scenario_init(&s, id);
    oracle_point_t pts[32];

    /* Print header for this scenario. */
    fprintf(stderr, "  --- %s ---\n", scenario_name(id));
    fprintf(stderr, "  frame |");
    for (int i = 0; i < n_ells; i++) fprintf(stderr, " %7.3f", ells[i]);
    fprintf(stderr, "    (ell)\n");

    int idx = 0;
    for (int t = 0; t < cfg.num_frames; t++) {
        sim_step(&s, NULL);
        if (idx < 3 && t == sample_frames[idx]) {
            oracle_eval_sweep(s.n, s.px, s.py, cfg.radius,
                              cfg.domain_w, cfg.domain_h,
                              ells, n_ells, 1, pts);
            fprintf(stderr, "  %5d |", t);
            for (int i = 0; i < n_ells; i++) fprintf(stderr, " %7ld", pts[i].S);
            fprintf(stderr, "    (S)\n");
            idx++;
        }
    }
    sim_free(&s);
}

static void print_S_tables_all(void) {
    fprintf(stderr, "\n--- (ell, S) sweep tables for power-law check ---\n");
    for (int id = 0; id < SCEN_COUNT; id++) {
        print_S_tables_one((scenario_id_t)id);
    }
}

/* =================================================================== main */

int main(void) {
    fprintf(stderr, "=== Phase 4 Verification Gates ===\n");

    int all_ok = 1;

    fprintf(stderr, "\n--- G1 oracle determinism ---\n");
    all_ok &= g1_oracle_deterministic();

    fprintf(stderr, "\n--- G2 non-perturbing ---\n");
    all_ok &= g2_non_perturbing();

    fprintf(stderr, "\n--- G3 sweep brackets the optimum (per scenario) ---\n");
    all_ok &= g3_sweep_brackets_all();

    fprintf(stderr, "\n--- G4 static oracle beats Ericson (per scenario) ---\n");
    all_ok &= g4_oracle_beats_ericson_all();

    print_S_tables_all();

    fprintf(stderr, "\n%s\n", all_ok ? "Phase 4: ALL OK" : "Phase 4: FAILED");
    return all_ok ? 0 : 1;
}
