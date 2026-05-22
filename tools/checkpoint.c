/*
 * Pre-Phase-5 data checkpoint.
 *
 * One-off scratch tool: uses the existing oracle + baselines APIs only,
 * no changes to src/. Not in Makefile; build manually:
 *
 *   gcc -O2 -march=native -Wall -Wextra -Wno-unused-parameter -std=c11 -Isrc \
 *       -o checkpoint tools/checkpoint.c \
 *       src/rng.c src/json.c src/config.c src/grid.c src/sim.c src/log.c \
 *       src/scenarios.c src/oracle.c src/baselines.c -lm
 *
 * Parts:
 *  2. Large-N (ell, S) spot-check on stable and multicluster at experiment
 *     scale -- to see if log S vs log ell is cleaner than at N=300.
 *  3. Per-frame vs static oracle headroom across all 8 scenarios.
 */

#define _POSIX_C_SOURCE 199309L

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

/* =================================================== Part 2: large-N spot-check */

static void largeN_spot_check(scenario_id_t id, int N, double W, double H,
                              int frames_to_run) {
    config_t cfg;
    config_default(&cfg);
    cfg.N           = N;
    cfg.domain_w    = W;  cfg.domain_h = H;
    cfg.radius      = 0.5;
    cfg.cell_size   = 2.0;       /* irrelevant for snapshot eval */
    cfg.dt          = 1.0 / 60.0;
    cfg.init_speed  = 5.0;
    cfg.num_frames  = frames_to_run;
    cfg.seed        = 42;
    cfg.restitution = 0.95;
    strncpy(cfg.scenario, scenario_name(id), sizeof(cfg.scenario) - 1);
    cfg.scenario[sizeof(cfg.scenario) - 1] = '\0';

    fprintf(stderr, "=== %s, N=%d on %.0fx%.0f ===\n",
            scenario_name(id), N, W, H);

    sim_t s; sim_init(&s, &cfg);
    fprintf(stderr, "  placing %d particles ...\n", N);
    scenario_init(&s, id);

    double ells[64];
    double ell_min = 2.0 * cfg.radius;
    double ell_max = 0.25 * fmin(W, H);
    int n_ells = oracle_make_geo_sweep(ell_min, ell_max, 1.25, ells, 64);

    fprintf(stderr, "  frame |");
    for (int i = 0; i < n_ells; i++) fprintf(stderr, " %9.2f", ells[i]);
    fprintf(stderr, "    (ell)\n");

    oracle_point_t pts[64];

    /* frame 1 */
    sim_step(&s, NULL);
    oracle_eval_sweep(s.n, s.px, s.py, cfg.radius, W, H, ells, n_ells, 1, pts);
    fprintf(stderr, "  %5d |", 1);
    for (int i = 0; i < n_ells; i++) fprintf(stderr, " %9ld", pts[i].S);
    fprintf(stderr, "    (S)\n");

    /* frame frames_to_run */
    for (int t = 1; t < frames_to_run; t++) sim_step(&s, NULL);
    oracle_eval_sweep(s.n, s.px, s.py, cfg.radius, W, H, ells, n_ells, 1, pts);
    fprintf(stderr, "  %5d |", frames_to_run);
    for (int i = 0; i < n_ells; i++) fprintf(stderr, " %9ld", pts[i].S);
    fprintf(stderr, "    (S)\n\n");

    sim_free(&s);
}

/* =================================================== Part 3: oracle headroom */

static void headroom_one(scenario_id_t id) {
    config_t cfg;
    config_default(&cfg);
    cfg.N           = 300;
    cfg.domain_w    = 60;  cfg.domain_h = 60;
    cfg.radius      = 0.5;
    cfg.dt          = 1.0 / 60.0;
    cfg.init_speed  = 10.0;
    cfg.num_frames  = 50;
    cfg.seed        = 31;
    cfg.restitution = 0.95;
    /* Configured cell_size for the trajectory: analytical baseline. */
    cfg.cell_size   = baseline_uniform_analytical(cfg.N, cfg.domain_w,
                                                  cfg.domain_h, cfg.radius);
    strncpy(cfg.scenario, scenario_name(id), sizeof(cfg.scenario) - 1);
    cfg.scenario[sizeof(cfg.scenario) - 1] = '\0';

    double ells[32];
    int n_ells = oracle_make_geo_sweep(2.0 * cfg.radius,
                                       0.5 * fmin(cfg.domain_w, cfg.domain_h),
                                       1.25, ells, 32);

    sim_t s; sim_init(&s, &cfg); scenario_init(&s, id);

    double **costs = malloc(sizeof(double *) * (size_t)cfg.num_frames);
    for (int t = 0; t < cfg.num_frames; t++)
        costs[t] = malloc(sizeof(double) * (size_t)n_ells);

    oracle_point_t pts[32];
    for (int t = 0; t < cfg.num_frames; t++) {
        sim_step(&s, NULL);
        oracle_eval_sweep(s.n, s.px, s.py, cfg.radius,
                          cfg.domain_w, cfg.domain_h,
                          ells, n_ells, 10, pts);   /* K=10 to suppress noise */
        for (int i = 0; i < n_ells; i++) costs[t][i] = pts[i].t_detect;
    }

    /* Per-frame oracle: min t_detect each frame. */
    double pf_sum = 0;
    for (int t = 0; t < cfg.num_frames; t++) {
        int best = 0;
        for (int i = 1; i < n_ells; i++)
            if (costs[t][i] < costs[t][best]) best = i;
        pf_sum += costs[t][best];
    }
    double pf_mean = pf_sum / (double)cfg.num_frames;

    /* Static (snapshot-based): best single ell summed across frames. */
    double best_sum = 1.0e18;
    int    best_idx = 0;
    for (int i = 0; i < n_ells; i++) {
        double sum = 0;
        for (int t = 0; t < cfg.num_frames; t++) sum += costs[t][i];
        if (sum < best_sum) { best_sum = sum; best_idx = i; }
    }
    double static_mean = best_sum / (double)cfg.num_frames;

    fprintf(stderr,
        "  %-12s | per-frame=%8.3f us | static(ell=%6.3f)=%8.3f us | "
        "pf/st=%.3f | headroom=%.0f%%\n",
        scenario_name(id),
        pf_mean * 1e6, ells[best_idx], static_mean * 1e6,
        pf_mean / static_mean,
        (1.0 - pf_mean / static_mean) * 100.0);

    for (int t = 0; t < cfg.num_frames; t++) free(costs[t]);
    free(costs);
    sim_free(&s);
}

/* ====================================================================== main */

int main(void) {
    fprintf(stderr, "==============================================================\n");
    fprintf(stderr, "Pre-Phase-5 data checkpoint  (no new functionality)\n");
    fprintf(stderr, "==============================================================\n\n");

    fprintf(stderr, "--- Part 2: large-N (ell, S) spot-check ---\n\n");
    largeN_spot_check(SCEN_STABLE, 50000, 800.0, 600.0, 5);
    fprintf(stderr,
        "  NOTE: multicluster's current cluster_R = min(W,H)*0.25*0.4 = 60\n"
        "  cannot pack 50000 particles at ~87%% per-cluster density (rejection\n"
        "  sampler would fail). Spot-checking at N=8000 (per-cluster ~14%%);\n"
        "  scaling cluster geometry with N is a follow-up before the experiment\n"
        "  campaign (already part of TODO.md item 1).\n\n");
    largeN_spot_check(SCEN_MULTICLUSTER, 8000, 800.0, 600.0, 5);

    fprintf(stderr,
        "--- Part 3: per-frame vs static oracle headroom ---\n\n"
        "  Both means computed from the SAME snapshot sweep (K=10 timing repeats).\n"
        "  Trajectory uses analytical-baseline cell_size. Per-frame oracle picks\n"
        "  best ell each frame; static (snapshot-based) picks single best ell over\n"
        "  all frames. Headroom = 1 - per_frame/static -- the upper bound on what\n"
        "  any adaptive rule could capture vs the best fixed cell size.\n\n");
    for (int id = 0; id < SCEN_COUNT; id++) {
        headroom_one((scenario_id_t)id);
    }

    return 0;
}
