/*
 * Pre-Phase-5 focused diagnostic.
 *
 * Long-run, large-N per-frame-vs-static oracle headroom + ell* range, on
 * dynamic scenarios (pulse, funnel, collapse, explosion) and two static
 * controls (stable, multicluster). Decides the paper's framing.
 *
 * One-off scratch tool, not in Makefile. Build:
 *   gcc -O2 -march=native -Wall -Wextra -Wno-unused-parameter -std=c11 -Isrc \
 *       -o diag tools/diag.c \
 *       src/rng.c src/json.c src/config.c src/grid.c src/sim.c src/log.c \
 *       src/scenarios.c src/oracle.c src/baselines.c -lm
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

typedef struct {
    scenario_id_t id;
    int           N;
    int           frames;

    /* outputs */
    double ell_min;
    double ell_max;
    double ell_ratio;
    double pf_mean_us;
    double static_ell;
    double static_mean_us;
    double headroom_pct;
    double pf_p95_us, pf_p99_us, pf_max_us;
    double st_p95_us, st_p99_us, st_max_us;
} result_t;

static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    return da < db ? -1 : (da > db ? 1 : 0);
}

/* Nearest-rank percentile on already-sorted array. */
static double pct_sorted(const double *sorted, int n, double q) {
    if (n <= 0) return 0;
    int k = (int)floor(q * (double)(n - 1) + 0.5);
    if (k < 0) k = 0;
    if (k >= n) k = n - 1;
    return sorted[k];
}

static void run_one(result_t *res) {
    config_t cfg;
    config_default(&cfg);
    cfg.N           = res->N;
    cfg.domain_w    = 800.0;
    cfg.domain_h    = 600.0;
    cfg.radius      = 0.5;
    cfg.dt          = 1.0 / 60.0;
    cfg.init_speed  = 10.0;
    cfg.num_frames  = res->frames;
    cfg.seed        = 42;
    cfg.restitution = 0.95;
    cfg.cell_size   = baseline_uniform_analytical(cfg.N, cfg.domain_w,
                                                  cfg.domain_h, cfg.radius);
    strncpy(cfg.scenario, scenario_name(res->id), sizeof(cfg.scenario) - 1);
    cfg.scenario[sizeof(cfg.scenario) - 1] = '\0';

    fprintf(stderr,
        "=== %s, N=%d, %d frames (config cell_size=%.2f, analytical) ===\n",
        scenario_name(res->id), cfg.N, cfg.num_frames, cfg.cell_size);
    fflush(stderr);

    sim_t s;
    sim_init(&s, &cfg);
    fprintf(stderr, "  placing %d particles ...\n", cfg.N);  fflush(stderr);
    scenario_init(&s, res->id);
    fprintf(stderr, "  running %d frames + oracle sweep ...\n", cfg.num_frames);
    fflush(stderr);

    /* Sweep -- ell in [2r, ~20] with gamma=1.25 (~14 points). Brackets the
     * expected optimum (sqrt(V/N) ≈ 3.1 at N=50K on 800x600) with margin
     * on both sides. ell_max=20 keeps per-eval pair tests bounded (avoids
     * the tens-of-millions-of-pairs regime that wastes runtime at large ell). */
    double ells[32];
    int n_ells = oracle_make_geo_sweep(2.0 * cfg.radius, 20.0, 1.25, ells, 32);

    double *t_per_ell    = (double *)malloc(sizeof(double) *
                              (size_t)cfg.num_frames * (size_t)n_ells);
    int    *best_ell_idx = (int    *)malloc(sizeof(int) * (size_t)cfg.num_frames);

    oracle_point_t pts[32];
    int progress_every = cfg.num_frames / 10;
    if (progress_every < 1) progress_every = 1;
    for (int t = 0; t < cfg.num_frames; t++) {
        sim_step(&s, NULL);
        oracle_eval_sweep(s.n, s.px, s.py, cfg.radius,
                          cfg.domain_w, cfg.domain_h,
                          ells, n_ells, 3, pts);
        int best = 0;
        for (int i = 1; i < n_ells; i++) {
            if (pts[i].t_detect < pts[best].t_detect) best = i;
        }
        best_ell_idx[t] = best;
        for (int i = 0; i < n_ells; i++) {
            t_per_ell[t * n_ells + i] = pts[i].t_detect;
        }
        if ((t + 1) % progress_every == 0) {
            fprintf(stderr, "    %d/%d frames  (last best ell=%.2f, t=%.1f us)\n",
                t + 1, cfg.num_frames, ells[best], pts[best].t_detect * 1e6);
            fflush(stderr);
        }
    }
    sim_free(&s);

    /* T1: per-frame oracle ell* min / max / ratio. */
    double e_min = ells[best_ell_idx[0]], e_max = e_min;
    for (int t = 1; t < cfg.num_frames; t++) {
        double e = ells[best_ell_idx[t]];
        if (e < e_min) e_min = e;
        if (e > e_max) e_max = e;
    }

    /* T2: mean detection cost. */
    double pf_sum = 0;
    for (int t = 0; t < cfg.num_frames; t++) {
        pf_sum += t_per_ell[t * n_ells + best_ell_idx[t]];
    }
    double pf_mean = pf_sum / (double)cfg.num_frames;

    double best_sum = 1e18; int best_idx = 0;
    for (int i = 0; i < n_ells; i++) {
        double sum = 0;
        for (int t = 0; t < cfg.num_frames; t++)
            sum += t_per_ell[t * n_ells + i];
        if (sum < best_sum) { best_sum = sum; best_idx = i; }
    }
    double static_mean = best_sum / (double)cfg.num_frames;

    /* T3: tail (p95, p99, max) of per-frame and static t_detect series. */
    double *pf_ser = (double *)malloc(sizeof(double) * (size_t)cfg.num_frames);
    double *st_ser = (double *)malloc(sizeof(double) * (size_t)cfg.num_frames);
    for (int t = 0; t < cfg.num_frames; t++) {
        pf_ser[t] = t_per_ell[t * n_ells + best_ell_idx[t]];
        st_ser[t] = t_per_ell[t * n_ells + best_idx];
    }
    qsort(pf_ser, cfg.num_frames, sizeof(double), cmp_double);
    qsort(st_ser, cfg.num_frames, sizeof(double), cmp_double);

    res->ell_min       = e_min;
    res->ell_max       = e_max;
    res->ell_ratio     = e_max / e_min;
    res->pf_mean_us    = pf_mean * 1e6;
    res->static_ell    = ells[best_idx];
    res->static_mean_us= static_mean * 1e6;
    res->headroom_pct  = (1.0 - pf_mean / static_mean) * 100.0;
    res->pf_p95_us = pct_sorted(pf_ser, cfg.num_frames, 0.95) * 1e6;
    res->pf_p99_us = pct_sorted(pf_ser, cfg.num_frames, 0.99) * 1e6;
    res->pf_max_us = pf_ser[cfg.num_frames - 1] * 1e6;
    res->st_p95_us = pct_sorted(st_ser, cfg.num_frames, 0.95) * 1e6;
    res->st_p99_us = pct_sorted(st_ser, cfg.num_frames, 0.99) * 1e6;
    res->st_max_us = st_ser[cfg.num_frames - 1] * 1e6;

    free(t_per_ell);
    free(best_ell_idx);
    free(pf_ser);
    free(st_ser);

    fprintf(stderr,
        "  RESULTS  ell* in [%.3f, %.3f] ratio=%.2f   pf=%.1f us  static(ell=%.2f)=%.1f us  headroom=%.0f%%\n\n",
        res->ell_min, res->ell_max, res->ell_ratio,
        res->pf_mean_us, res->static_ell, res->static_mean_us, res->headroom_pct);
}

int main(void) {
    fprintf(stderr,
        "===================================================================\n"
        "Pre-Phase-5 diagnostic: per-frame vs static oracle (large N, long runs)\n"
        "===================================================================\n\n");

    result_t res[] = {
        { SCEN_PULSE,        50000, 600, 0,0,0,0,0,0,0,0,0,0,0,0,0 },
        { SCEN_FUNNEL,       50000, 600, 0,0,0,0,0,0,0,0,0,0,0,0,0 },
        { SCEN_COLLAPSE,     50000, 400, 0,0,0,0,0,0,0,0,0,0,0,0,0 },
        { SCEN_EXPLOSION,    50000, 400, 0,0,0,0,0,0,0,0,0,0,0,0,0 },
        { SCEN_STABLE,       50000, 400, 0,0,0,0,0,0,0,0,0,0,0,0,0 },
        /* multicluster: current cluster_R cannot fit N=50K (~87% per-cluster
         * packing); the largest that fills cleanly at the current geometry is
         * ~20K (per-cluster ~14%). Recorded in TODO.md item 1.            */
        { SCEN_MULTICLUSTER, 20000, 400, 0,0,0,0,0,0,0,0,0,0,0,0,0 },
    };
    int n = (int)(sizeof(res) / sizeof(res[0]));

    for (int i = 0; i < n; i++) run_one(&res[i]);

    /* ============== Summary tables ============== */
    fprintf(stderr, "===================================================================\n");
    fprintf(stderr, "Summary\n");
    fprintf(stderr, "===================================================================\n\n");

    fprintf(stderr, "Table 1: per-frame oracle ell* range\n");
    fprintf(stderr, "  scenario     |    N   | frames |  ell_min  |  ell_max  |  ratio\n");
    fprintf(stderr, "  -------------+--------+--------+-----------+-----------+--------\n");
    for (int i = 0; i < n; i++) {
        fprintf(stderr, "  %-12s | %6d | %6d |  %8.3f |  %8.3f |  %5.2f\n",
            scenario_name(res[i].id), res[i].N, res[i].frames,
            res[i].ell_min, res[i].ell_max, res[i].ell_ratio);
    }

    fprintf(stderr,
        "\nTable 2: mean headroom (per-frame oracle vs static-best fixed ell)\n");
    fprintf(stderr,
        "  scenario     |   pf mean   | static ell |  static mean |  pf/static |  headroom\n");
    fprintf(stderr,
        "  -------------+-------------+------------+--------------+------------+------------\n");
    for (int i = 0; i < n; i++) {
        fprintf(stderr,
            "  %-12s |  %8.1f us |  %8.3f  |   %8.1f us |   %6.3f   |   %5.1f%%\n",
            scenario_name(res[i].id),
            res[i].pf_mean_us, res[i].static_ell, res[i].static_mean_us,
            res[i].pf_mean_us / res[i].static_mean_us, res[i].headroom_pct);
    }

    fprintf(stderr, "\nTable 3: tail headroom (us; ratio = pf / static)\n");
    fprintf(stderr,
        "  scenario     |   p95 pf   |   p95 st   |   ratio  |   p99 pf   |   p99 st   |   ratio  |   max pf   |   max st   |   ratio\n");
    fprintf(stderr,
        "  -------------+------------+------------+----------+------------+------------+----------+------------+------------+---------\n");
    for (int i = 0; i < n; i++) {
        fprintf(stderr,
            "  %-12s |  %8.1f  |  %8.1f  |  %6.3f  |  %8.1f  |  %8.1f  |  %6.3f  |  %8.1f  |  %8.1f  |  %6.3f\n",
            scenario_name(res[i].id),
            res[i].pf_p95_us, res[i].st_p95_us, res[i].pf_p95_us / res[i].st_p95_us,
            res[i].pf_p99_us, res[i].st_p99_us, res[i].pf_p99_us / res[i].st_p99_us,
            res[i].pf_max_us, res[i].st_max_us, res[i].pf_max_us / res[i].st_max_us);
    }

    return 0;
}
