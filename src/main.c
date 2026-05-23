#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "config.h"
#include "sim.h"
#include "scenarios.h"
#include "log.h"
#include "oracle.h"
#include "adapter.h"
#include <string.h>

#define ORACLE_MAX_PTS 64

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void dump_snapshot(FILE *snap, int frame, const sim_t *s) {
    if (!snap) return;
    for (int i = 0; i < s->n; i++) {
        fprintf(snap, "%d,%d,%.17g,%.17g,%.17g,%.17g\n",
                frame, i, s->px[i], s->py[i], s->vx[i], s->vy[i]);
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s <config.json>\n", argv[0]);
        return 1;
    }
    config_t cfg;
    if (config_load(argv[1], &cfg) != 0) return 1;
    config_print(&cfg);

    sim_t s;
    sim_init(&s, &cfg);
    scenario_id_t scen = scenario_from_name(cfg.scenario);
    scenario_init(&s, scen);
    fprintf(stderr,
        "scenario: %s (forces=%s)\n",
        scenario_name(scen), scenario_has_forces(scen) ? "yes" : "no");

    FILE *csv = NULL;
    if (cfg.log_enabled && cfg.log_csv_path[0] != '\0') {
        csv = log_open(cfg.log_csv_path);
        if (csv) log_write_header(csv);
    }

    /* Optional per-frame oracle (Phase 4): non-perturbing sweep over a
     * geometric ladder of candidate cell sizes, evaluated on the frozen
     * post-step snapshot. Writes one row per (frame, ell) pair. */
    FILE  *oracle_csv = NULL;
    double oracle_ells[ORACLE_MAX_PTS];
    int    n_oracle_ells = 0;
    int    oracle_K_use  = 5;
    if (cfg.oracle_csv_path[0] != '\0') {
        double ell_min = cfg.oracle_ell_min > 0.0 ? cfg.oracle_ell_min : 2.0 * cfg.radius;
        double ell_max = cfg.oracle_ell_max > 0.0 ? cfg.oracle_ell_max
                                                 : 0.5 * fmin(cfg.domain_w, cfg.domain_h);
        double gamma   = cfg.oracle_gamma   > 1.0 ? cfg.oracle_gamma   : 1.25;
        if (cfg.oracle_K > 0) oracle_K_use = cfg.oracle_K;
        n_oracle_ells = oracle_make_geo_sweep(ell_min, ell_max, gamma,
                                              oracle_ells, ORACLE_MAX_PTS);
        oracle_csv = fopen(cfg.oracle_csv_path, "w");
        if (oracle_csv) {
            fputs("frame,ell,num_cells,S,candidate_pairs,"
                  "t_grid,t_broad,t_narrow,t_detect\n", oracle_csv);
            fprintf(stderr,
                "oracle: %d candidate ells in [%.3f, %.3f], gamma=%.3f, K=%d -> %s\n",
                n_oracle_ells, ell_min, ell_max, gamma, oracle_K_use,
                cfg.oracle_csv_path);
        }
    }

    /* Optional position snapshots: 5 frames evenly spaced. */
    FILE *snap = NULL;
    int snap_frames[5] = {-1, -1, -1, -1, -1};
    if (cfg.snapshot_csv_path[0] != '\0') {
        snap = fopen(cfg.snapshot_csv_path, "w");
        if (snap) {
            fputs("frame,id,x,y,vx,vy\n", snap);
            for (int k = 0; k < 5; k++) {
                int f = (cfg.num_frames * k) / 5;
                if (f >= cfg.num_frames) f = cfg.num_frames - 1;
                snap_frames[k] = f;
            }
        } else {
            fprintf(stderr, "main: failed to open snapshot '%s'\n", cfg.snapshot_csv_path);
        }
    }

    /* Phase 5: optional adaptive cell-size controller. */
    adapter_t     adapter;
    adapter_mode_t adapter_mode = ADAPTER_MODE_OFF;
    FILE         *adapter_csv   = NULL;
    int           adapter_on    = 0;
    if (cfg.adapter_mode[0] != '\0') {
        if      (strcmp(cfg.adapter_mode, "modeA") == 0) adapter_mode = ADAPTER_MODE_A;
        else if (strcmp(cfg.adapter_mode, "modeB") == 0) adapter_mode = ADAPTER_MODE_B;
        else if (strcmp(cfg.adapter_mode, "blind") == 0) adapter_mode = ADAPTER_MODE_BLIND;
        else if (strcmp(cfg.adapter_mode, "off")   == 0) adapter_mode = ADAPTER_MODE_OFF;
        else fprintf(stderr, "main: unknown adapter_mode '%s' -- defaulting to off\n",
                     cfg.adapter_mode);
        if (adapter_mode != ADAPTER_MODE_OFF) {
            adapter_init_default(&adapter, adapter_mode, 2,
                                 cfg.domain_w, cfg.domain_h,
                                 cfg.cell_size, cfg.radius);
            if (cfg.adapter_K_D2    > 0)   adapter.K_D2     = cfg.adapter_K_D2;
            if (cfg.adapter_gamma   > 1.0) adapter.gamma    = cfg.adapter_gamma;
            if (cfg.adapter_alpha_w > 0.0) adapter.alpha_w  = cfg.adapter_alpha_w;
            if (cfg.adapter_alpha_D > 0.0) adapter.alpha_D  = cfg.adapter_alpha_D;
            if (cfg.adapter_eps     > 0.0) adapter.eps      = cfg.adapter_eps;
            if (cfg.adapter_ell_max > 0.0) adapter.ell_max  = cfg.adapter_ell_max;
            adapter_on = 1;
            fprintf(stderr,
                "adapter: mode=%s ell0=%.3f ell_min=%.3f ell_max=%.3f "
                "K_D2=%d gamma=%.2f alpha_w=%.2f alpha_D=%.2f eps=%.3f\n",
                cfg.adapter_mode, adapter.ell_cur, adapter.ell_min, adapter.ell_max,
                adapter.K_D2, adapter.gamma, adapter.alpha_w, adapter.alpha_D, adapter.eps);
        }
    }
    if (adapter_on && cfg.adapter_csv_path[0] != '\0') {
        adapter_csv = fopen(cfg.adapter_csv_path, "w");
        if (adapter_csv) fputs(
            "frame,ell_cur,ell_raw,D2_hat,a_hat,b_hat,M_cur,S_cur,t_M,t_S,"
            "resized,did_probe,overhead_s\n", adapter_csv);
    }

    double ke0 = sim_total_kinetic_energy(&s);
    double pxt0, pyt0;
    sim_total_momentum(&s, &pxt0, &pyt0);

    double t0 = now_seconds();
    long total_pairs = 0, total_cands = 0;
    oracle_point_t oracle_pts[ORACLE_MAX_PTS];
    for (int t = 0; t < cfg.num_frames; t++) {
        sim_metrics_t m = { .frame = t };
        /* Always pass &m -- the adapter needs S_cur / M_cur / t_M / t_S
         * whether or not CSV logging is enabled. */
        sim_step(&s, &m);
        if (csv) log_write_row(csv, &m);

        if (adapter_on) {
            /* Spec §7: t_M = clear + prefix + cell-iteration (broad-phase
             * outer loop); t_S = pair-test inner work + narrow phase. */
            double t_S_in = (m.t_broad - m.t_broad_iter) + m.t_narrow;
            if (t_S_in < 0) t_S_in = 0;            /* timer noise floor   */
            double t_M_in = m.t_M + m.t_broad_iter;
            double ell_next;
            if (adapter_mode == ADAPTER_MODE_BLIND) {
                ell_next = adapter_blind_step(&adapter, s.n, s.px, s.py);
            } else {
                ell_next = adapter_step(&adapter, t_M_in, t_S_in, m.S, m.num_cells,
                                        s.n, s.px, s.py);
            }
            if (adapter_csv) fprintf(adapter_csv,
                "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%ld,%.9g,%.9g,%d,%d,%.9g\n",
                t, adapter.ell_cur, adapter.last_ell_raw, adapter.D2_hat,
                adapter.a_hat, adapter.b_hat, m.num_cells, m.S,
                t_M_in, t_S_in,
                adapter.last_resized, adapter.last_did_probe,
                adapter.last_overhead_s);
            if (adapter.last_resized) grid_resize(&s.grid, ell_next);
        }
        for (int k = 0; k < 5; k++) {
            if (t == snap_frames[k]) { dump_snapshot(snap, t, &s); break; }
        }
        /* Per-frame oracle, on the frozen post-step snapshot.
         * MUST be after sim_step (so positions reflect frame t) and MUST
         * not mutate s -- oracle_eval allocates its own grid/pair buffers. */
        if (oracle_csv && n_oracle_ells > 0) {
            oracle_eval_sweep(s.n, s.px, s.py, cfg.radius,
                              cfg.domain_w, cfg.domain_h,
                              oracle_ells, n_oracle_ells, oracle_K_use,
                              oracle_pts);
            for (int i = 0; i < n_oracle_ells; i++) {
                fprintf(oracle_csv,
                        "%d,%.17g,%d,%ld,%d,%.9g,%.9g,%.9g,%.9g\n",
                        t, oracle_pts[i].ell, oracle_pts[i].num_cells,
                        oracle_pts[i].S, oracle_pts[i].candidate_pairs,
                        oracle_pts[i].t_grid, oracle_pts[i].t_broad,
                        oracle_pts[i].t_narrow, oracle_pts[i].t_detect);
            }
        }
        total_pairs += s.last_pair_count;
        total_cands += s.last_cand_count;
    }
    double t1 = now_seconds();

    log_close(csv);
    if (snap) fclose(snap);
    if (oracle_csv) fclose(oracle_csv);
    if (adapter_csv) fclose(adapter_csv);

    double ke1 = sim_total_kinetic_energy(&s);
    double pxt1, pyt1;
    sim_total_momentum(&s, &pxt1, &pyt1);

    fprintf(stderr,
        "frames=%d  wall=%.3fs  ms/frame=%.4f  candidates=%ld  collisions=%ld  M=%d  S=%ld\n",
        cfg.num_frames, t1 - t0, (t1 - t0) * 1000.0 / cfg.num_frames,
        total_cands, total_pairs, grid_M(&s.grid), grid_S(&s.grid));
    fprintf(stderr,
        "KE: initial=%.6g final=%.6g  ratio=%.9f\n", ke0, ke1, ke1 / ke0);
    fprintf(stderr,
        "P : initial=(%.6g, %.6g)  final=(%.6g, %.6g)\n", pxt0, pyt0, pxt1, pyt1);
    if (csv != NULL) fprintf(stderr, "CSV log:        %s\n", cfg.log_csv_path);
    if (snap_frames[0] >= 0 && cfg.snapshot_csv_path[0])
        fprintf(stderr, "snapshot CSV:   %s (frames %d,%d,%d,%d,%d)\n",
                cfg.snapshot_csv_path,
                snap_frames[0], snap_frames[1], snap_frames[2],
                snap_frames[3], snap_frames[4]);

    sim_free(&s);
    return 0;
}
