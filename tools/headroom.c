/*
 * Self-validated headroom measurement (Phase 7, Task 1).
 *
 * Headroom = how much a per-frame-optimal ell beats the best single fixed ell.
 * The earlier compare harness measured the two oracles with DIFFERENT primitives
 * (per-frame oracle via oracle_eval on fresh grids; static oracle via the sim's
 * own sim_step path, which is single-shot and even folds the broad-phase
 * iter-time dryrun into t_narrow), so even a genuinely static distribution read
 * a spurious ~13% headroom.
 *
 * Here BOTH oracles are derived from ONE shared (frame × ell) cost matrix,
 * measured by the SAME primitive (oracle_eval_split, de-biased: select on one
 * timing sample, report the other):
 *   - per-frame oracle = mean_f  eval[f][ argmin_e sel[f][e] ]
 *   - static  oracle   = mean_f  eval[f][ e* ],  e* = argmin_e mean_f sel[f][e]
 * They cannot diverge from a measurement inconsistency, only from genuine
 * distribution change across frames.
 *
 * FROZEN control: feed the identical mid-run snapshot as every frame. Then
 * every row of the matrix is the same distribution, so per-frame and static
 * MUST coincide -- the frozen headroom must be ~0 within timing noise. That is
 * the measurement's self-test.
 *
 * Usage: ./headroom <N> [frames]   (runs frozen controls + live, 8 2D scenarios)
 */

#define _POSIX_C_SOURCE 199309L
#include "config.h"
#include "sim.h"
#include "scenarios.h"
#include "grid.h"
#include "oracle.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double pf_mean_us;      /* per-frame oracle cost (mean over frames)   */
    double static_mean_us;  /* static oracle cost (best single ell)       */
    double static_ell;
    double pf_ell_mean;
    double headroom_pct;    /* (static/pf - 1) * 100                       */
    int    n_frames;
} headroom_t;

static config_t base_cfg(scenario_id_t id, int N) {
    config_t cfg; config_default(&cfg);
    cfg.dim = 2;
    cfg.N = N;
    cfg.domain_w = 800.0; cfg.domain_h = 600.0;
    cfg.radius = 0.5; cfg.dt = 1.0/60.0;
    cfg.init_speed = 10.0; cfg.restitution = 0.95;
    cfg.num_frames = 300; cfg.seed = 1;
    cfg.cell_size = 2.0;
    snprintf(cfg.scenario, sizeof(cfg.scenario), "%s", scenario_name(id));
    return cfg;
}

/* Build the shared (F × n_ells) sel/eval matrices and derive both oracles.
 * frozen!=0 -> capture one mid-run snapshot and reuse it for every frame. */
static headroom_t measure(scenario_id_t id, int N, int frozen,
                          const double *ells, int n_ells, int K, int F) {
    config_t cfg = base_cfg(id, N);
    cfg.log_enabled = 0; cfg.log_csv_path[0]='\0';
    cfg.snapshot_csv_path[0]='\0'; cfg.oracle_csv_path[0]='\0'; cfg.adapter_mode[0]='\0';

    sim_t s; sim_init(&s, &cfg);
    scenario_init(&s, id);

    double *fx = NULL, *fy = NULL;
    if (frozen) {
        for (int t = 0; t < F / 2; t++) sim_step(&s, NULL);   /* mid-run config */
        fx = (double *)malloc(sizeof(double) * s.n);
        fy = (double *)malloc(sizeof(double) * s.n);
        memcpy(fx, s.px, sizeof(double) * s.n);
        memcpy(fy, s.py, sizeof(double) * s.n);
    }

    double *sel  = (double *)malloc(sizeof(double) * (size_t)F * n_ells);
    double *eval = (double *)malloc(sizeof(double) * (size_t)F * n_ells);

    for (int f = 0; f < F; f++) {
        const double *px, *py;
        if (frozen) { px = fx; py = fy; }            /* identical every frame */
        else        { sim_step(&s, NULL); px = s.px; py = s.py; }
        for (int e = 0; e < n_ells; e++) {
            oracle_point_t os, oe;
            oracle_eval_split(s.n, px, py, cfg.radius, cfg.domain_w, cfg.domain_h,
                              ells[e], K, &os, &oe);
            sel [f * n_ells + e] = os.t_detect;
            eval[f * n_ells + e] = oe.t_detect;
        }
    }

    /* per-frame oracle: select per-frame on sel, report eval */
    double pf = 0, pf_ell = 0;
    for (int f = 0; f < F; f++) {
        int b = 0;
        for (int e = 1; e < n_ells; e++)
            if (sel[f * n_ells + e] < sel[f * n_ells + b]) b = e;
        pf     += eval[f * n_ells + b];
        pf_ell += ells[b];
    }
    pf /= F; pf_ell /= F;

    /* static oracle: select single best ell on mean sel, report mean eval */
    int best_e = 0; double best_mean = 1e300;
    for (int e = 0; e < n_ells; e++) {
        double m = 0;
        for (int f = 0; f < F; f++) m += sel[f * n_ells + e];
        m /= F;
        if (m < best_mean) { best_mean = m; best_e = e; }
    }
    double st = 0;
    for (int f = 0; f < F; f++) st += eval[f * n_ells + best_e];
    st /= F;

    headroom_t h;
    h.pf_mean_us = pf * 1e6;
    h.static_mean_us = st * 1e6;
    h.static_ell = ells[best_e];
    h.pf_ell_mean = pf_ell;
    h.headroom_pct = (st / pf - 1.0) * 100.0;
    h.n_frames = F;

    free(sel); free(eval); free(fx); free(fy);
    sim_free(&s);
    return h;
}

int main(int argc, char **argv) {
    int N = 8000, F = 40, K = 2;
    if (argc > 1) N = atoi(argv[1]);
    if (argc > 2) F = atoi(argv[2]);

    /* clipped operative ladder [2r, 20] */
    double ells[32];
    int n_ells = oracle_make_geo_sweep(1.0, 20.0, 1.3, ells, 32);

    fprintf(stderr, "=== Headroom (N=%d, F=%d, K=%d/group, ladder [%.1f,%.1f] %d pts) ===\n",
            N, F, K, ells[0], ells[n_ells-1], n_ells);

    /* ---- FROZEN CONTROL (self-test): must read ~0 ---- */
    fprintf(stderr, "\n-- frozen control (identical snapshot every frame; headroom must be ~0) --\n");
    scenario_id_t frozen_ids[] = { SCEN_STABLE, SCEN_VORTEX };
    for (int i = 0; i < 2; i++) {
        headroom_t h = measure(frozen_ids[i], N, 1, ells, n_ells, K, F);
        fprintf(stderr, "  %-12s : pf=%.1f us  static=%.1f us  static_ell=%.2f  headroom=%+.2f%%\n",
                scenario_name(frozen_ids[i]), h.pf_mean_us, h.static_mean_us, h.static_ell, h.headroom_pct);
    }

    /* ---- LIVE headroom, all 8 2D scenarios ---- */
    fprintf(stderr, "\n-- live headroom (per-frame oracle vs best fixed ell) --\n");
    fprintf(stderr, "  %-12s |    pf_us | static_us | pf_ell | stat_ell | headroom\n", "scenario");
    fprintf(stderr, "  -------------+----------+-----------+--------+----------+---------\n");
    for (int id = 0; id < SCEN_COUNT; id++) {
        if (id == SCEN_SETTLING) continue;            /* 8 spec 2D scenarios */
        if (id == SCEN_MULTICLUSTER && N > 25000) {
            fprintf(stderr, "  %-12s |  (skipped: placement saturates at N=%d in 800x600)\n",
                    scenario_name((scenario_id_t)id), N);
            continue;
        }
        headroom_t h = measure((scenario_id_t)id, N, 0, ells, n_ells, K, F);
        fprintf(stderr, "  %-12s | %8.1f | %9.1f | %6.2f | %8.2f | %+7.2f%%\n",
                scenario_name((scenario_id_t)id),
                h.pf_mean_us, h.static_mean_us, h.pf_ell_mean, h.static_ell, h.headroom_pct);
    }
    return 0;
}
