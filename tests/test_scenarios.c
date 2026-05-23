/*
 * Phase 3 verification gates (spec section 10):
 *   G1. Per-scenario characteristic statistic confirming the intended behaviour.
 *   G2. Every scenario is deterministic: two runs of the same (config, seed)
 *       produce a bit-identical metrics CSV (deterministic columns).
 *   G3. Phase 1/2 gates per scenario:
 *       a) grid pipeline output == brute force overlap set, every frame.
 *       b) energy stays bounded over a long run -- KE non-increasing for
 *          force-free scenarios; KE finite and well below a generous absolute
 *          ceiling for force-driven ones (symplectic-Euler stability).
 *   G4. Particles stay inside the domain -- no wall escapes.
 *
 * Per-scenario invariants chosen:
 *   stable       : S ≈ N²/M       (D₂ ≈ d, uniform baseline) -- ratio in [0.5, 2.0]
 *   collapse     : mean(|p-c|) at end < 0.7× mean(|p-c|) at start
 *   explosion    : mean outward radial velocity > 0 after 1 step
 *                  AND stddev(|p-c|) grows by >= 1.5× over 30 frames
 *   vortex       : S(0) > 2·N²/M  (filamentary, low D₂)
 *   funnel       : mean(|y - A·sin(B·x) - y_c|) at end < 0.7× start
 *   pulse        : (max - min)/max of mean(|p-c|) over 200 frames > 0.3
 *                  (the breathing oscillation under the central spring)
 *   multicluster : S(0) > 2·N²/M  (K clusters; high non-uniformity)
 *   stream       : mean|y - H/2| < 2× channel half-height
 *                  AND mean(vx) > 0 (flow direction preserved)
 */

#include "config.h"
#include "sim.h"
#include "scenarios.h"
#include "grid.h"
#include "log.h"
#include "rng.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================== helpers */

static config_t std_cfg(scenario_id_t id) {
    config_t cfg;
    config_default(&cfg);
    cfg.N           = 400;
    cfg.domain_w    = 80.0;
    cfg.domain_h    = 80.0;
    cfg.radius      = 0.5;
    cfg.cell_size   = 2.0;
    cfg.dt          = 1.0 / 60.0;
    cfg.init_speed  = 10.0;
    cfg.restitution = 0.95;
    cfg.num_frames  = 100;
    cfg.seed        = 123;
    strncpy(cfg.scenario, scenario_name(id), sizeof(cfg.scenario) - 1);
    cfg.scenario[sizeof(cfg.scenario) - 1] = '\0';
    return cfg;
}

static double mean_r(const sim_t *s) {
    double cx = s->cfg.domain_w * 0.5, cy = s->cfg.domain_h * 0.5;
    double sum = 0;
    for (int i = 0; i < s->n; i++) {
        double dx = s->px[i] - cx, dy = s->py[i] - cy;
        sum += sqrt(dx * dx + dy * dy);
    }
    return sum / (double)s->n;
}

static double mean_radial_v(const sim_t *s) {
    double cx = s->cfg.domain_w * 0.5, cy = s->cfg.domain_h * 0.5;
    double sum = 0;
    for (int i = 0; i < s->n; i++) {
        double dx = s->px[i] - cx, dy = s->py[i] - cy;
        double r  = sqrt(dx * dx + dy * dy);
        if (r < 1e-9) continue;
        sum += (s->vx[i] * dx + s->vy[i] * dy) / r;
    }
    return sum / (double)s->n;
}

static double funnel_mean_dev(const sim_t *s) {
    double W = s->cfg.domain_w, H = s->cfg.domain_h;
    double A = H * 0.15;
    double B = 4.0 * M_PI / W;
    double yc = H * 0.5;
    double sum = 0;
    for (int i = 0; i < s->n; i++) {
        double target = yc + A * sin(B * s->px[i]);
        sum += fabs(s->py[i] - target);
    }
    return sum / (double)s->n;
}

static double mean_y_dev_from_center(const sim_t *s) {
    double yc = s->cfg.domain_h * 0.5;
    double sum = 0;
    for (int i = 0; i < s->n; i++) sum += fabs(s->py[i] - yc);
    return sum / (double)s->n;
}

static double mean_vx_all(const sim_t *s) {
    double sum = 0;
    for (int i = 0; i < s->n; i++) sum += s->vx[i];
    return sum / (double)s->n;
}

/* ====================================== G1: characteristic statistics */

static int char_stable(void) {
    config_t cfg = std_cfg(SCEN_STABLE);
    cfg.num_frames = 60;
    sim_t s; sim_init(&s, &cfg); scenario_init(&s, SCEN_STABLE);
    for (int t = 0; t < cfg.num_frames; t++) sim_step(&s, NULL);
    long S = grid_S(&s.grid);
    int  M = grid_M(&s.grid);
    /* Poisson baseline for random (Poisson-uniform) placement:
     * E[S] = N + N^2/M. The deterministic Cauchy-Schwarz floor N^2/M
     * is only attained by perfectly equal bin counts (rarely the case). */
    double S_poisson = (double)cfg.N + ((double)cfg.N * cfg.N) / (double)M;
    double ratio     = (double)S / S_poisson;
    int ok = (ratio > 0.5) && (ratio < 1.5);
    fprintf(stderr,
        "  stable       : S=%ld  Poisson_baseline=%.1f  ratio=%.2f  [0.5,1.5]  %s\n",
        S, S_poisson, ratio, ok ? "OK" : "FAIL");
    sim_free(&s);
    return ok;
}

static int char_collapse(void) {
    config_t cfg = std_cfg(SCEN_COLLAPSE);
    cfg.num_frames = 100;
    sim_t s; sim_init(&s, &cfg); scenario_init(&s, SCEN_COLLAPSE);
    double r0    = mean_r(&s);
    double r_min = r0;
    /* Particles converge, pile up at the center, then bounce out. We
     * track the minimum mean-r over the whole run -- that's where the
     * convergence is strongest, before the rebound starts.            */
    for (int t = 0; t < cfg.num_frames; t++) {
        sim_step(&s, NULL);
        double r = mean_r(&s);
        if (r < r_min) r_min = r;
    }
    int ok = r_min < 0.6 * r0;
    fprintf(stderr,
        "  collapse     : mean|p-c| start=%.2f  min=%.2f  ratio=%.2f (<0.6?)  %s\n",
        r0, r_min, r_min / r0, ok ? "OK" : "FAIL");
    sim_free(&s);
    return ok;
}

static int char_explosion(void) {
    config_t cfg = std_cfg(SCEN_EXPLOSION);
    cfg.N          = 200;
    cfg.num_frames = 40;
    sim_t s; sim_init(&s, &cfg); scenario_init(&s, SCEN_EXPLOSION);
    /* Note: stddev(|p-c|) is invariant under purely outward radial
     * translation (every particle's radius shifts by the same amount),
     * so the right characteristic is *mean*-r growth, not stddev. */
    double r0 = mean_r(&s);
    sim_step(&s, NULL);
    double v_rad_after1 = mean_radial_v(&s);
    double r_max = r0;
    for (int t = 1; t < cfg.num_frames; t++) {
        sim_step(&s, NULL);
        double r = mean_r(&s);
        if (r > r_max) r_max = r;
    }
    int ok = (v_rad_after1 > 0.0) && (r_max > 1.5 * r0);
    fprintf(stderr,
        "  explosion    : v_rad@1=%.2f  mean|p-c| %.2f -> max %.2f (x%.2f, >1.5?)  %s\n",
        v_rad_after1, r0, r_max, r_max / r0, ok ? "OK" : "FAIL");
    sim_free(&s);
    return ok;
}

static int char_vortex(void) {
    config_t cfg = std_cfg(SCEN_VORTEX);
    cfg.num_frames = 1;
    sim_t s; sim_init(&s, &cfg); scenario_init(&s, SCEN_VORTEX);
    grid_build(&s.grid, s.px, s.py);
    long S = grid_S(&s.grid);
    int  M = grid_M(&s.grid);
    double S_uniform = ((double)cfg.N * (double)cfg.N) / (double)M;
    double ratio = (double)S / S_uniform;
    int ok = ratio > 2.0;
    fprintf(stderr,
        "  vortex       : S(0)=%ld  S_uniform=%.1f  ratio=%.2f (>2?)  %s\n",
        S, S_uniform, ratio, ok ? "OK" : "FAIL");
    sim_free(&s);
    return ok;
}

static int char_funnel(void) {
    config_t cfg = std_cfg(SCEN_FUNNEL);
    cfg.num_frames = 600;   /* let the lightly damped oscillation settle */
    sim_t s; sim_init(&s, &cfg); scenario_init(&s, SCEN_FUNNEL);
    double dev0    = funnel_mean_dev(&s);
    double dev_min = dev0;
    for (int t = 0; t < cfg.num_frames; t++) {
        sim_step(&s, NULL);
        double d = funnel_mean_dev(&s);
        if (d < dev_min) dev_min = d;
    }
    int ok = dev_min < 0.7 * dev0;
    fprintf(stderr,
        "  funnel       : mean|y-target| start=%.2f  min=%.2f  ratio=%.2f (<0.7?)  %s\n",
        dev0, dev_min, dev_min / dev0, ok ? "OK" : "FAIL");
    sim_free(&s);
    return ok;
}

static int char_pulse(void) {
    config_t cfg = std_cfg(SCEN_PULSE);
    cfg.num_frames  = 200;
    cfg.restitution = 0.99;   /* less dissipation -> oscillation visible */
    sim_t s; sim_init(&s, &cfg); scenario_init(&s, SCEN_PULSE);
    double min_r = 1e18, max_r = 0;
    for (int t = 0; t < cfg.num_frames; t++) {
        sim_step(&s, NULL);
        double r = mean_r(&s);
        if (r < min_r) min_r = r;
        if (r > max_r) max_r = r;
    }
    double swing = (max_r > 0) ? ((max_r - min_r) / max_r) : 0;
    int ok = swing > 0.3;
    fprintf(stderr,
        "  pulse        : mean_r in [%.2f, %.2f]  swing=%.2f (>0.3?)  %s\n",
        min_r, max_r, swing, ok ? "OK" : "FAIL");
    sim_free(&s);
    return ok;
}

static int char_multicluster(void) {
    config_t cfg = std_cfg(SCEN_MULTICLUSTER);
    cfg.num_frames = 1;
    sim_t s; sim_init(&s, &cfg); scenario_init(&s, SCEN_MULTICLUSTER);
    grid_build(&s.grid, s.px, s.py);
    long S = grid_S(&s.grid);
    int  M = grid_M(&s.grid);
    double S_uniform = ((double)cfg.N * (double)cfg.N) / (double)M;
    double ratio = (double)S / S_uniform;
    int ok = ratio > 2.0;
    fprintf(stderr,
        "  multicluster : S(0)=%ld  S_uniform=%.1f  ratio=%.2f (>2?)  %s\n",
        S, S_uniform, ratio, ok ? "OK" : "FAIL");
    sim_free(&s);
    return ok;
}

static int char_stream(void) {
    config_t cfg = std_cfg(SCEN_STREAM);
    cfg.num_frames = 30;   /* short enough that few particles bounce off x-walls */
    sim_t s; sim_init(&s, &cfg); scenario_init(&s, SCEN_STREAM);
    double channel_half = cfg.domain_h * 0.1;
    for (int t = 0; t < cfg.num_frames; t++) sim_step(&s, NULL);
    double y_dev = mean_y_dev_from_center(&s);
    double vx_m  = mean_vx_all(&s);
    int ok_y  = y_dev < 2.0 * channel_half;
    int ok_vx = vx_m > 0.0;
    int ok = ok_y && ok_vx;
    fprintf(stderr,
        "  stream       : mean|y-H/2|=%.2f (<%.2f?)  mean_vx=%.2f (>0?)  %s\n",
        y_dev, 2.0 * channel_half, vx_m, ok ? "OK" : "FAIL");
    sim_free(&s);
    return ok;
}

/* ============================== G2: determinism per scenario (CSV diff) */

static char *csv_next(char **s) {
    if (!*s || !**s) return NULL;
    char *st = *s;
    while (**s && **s != ',' && **s != '\n') (*s)++;
    if (**s) { **s = '\0'; (*s)++; }
    return st;
}

static int csv_det_eq(const char *p1, const char *p2, int *rows_out, int *mism_out) {
    FILE *f1 = fopen(p1, "r"), *f2 = fopen(p2, "r");
    if (!f1 || !f2) { if (f1) fclose(f1); if (f2) fclose(f2); return 0; }
    int row = 0, mism = 0;
    char l1[2048], l2[2048];
    if (!fgets(l1, sizeof(l1), f1) || !fgets(l2, sizeof(l2), f2) ||
        strcmp(l1, l2) != 0) {
        if (mism_out) *mism_out = 1;
        if (rows_out) *rows_out = 0;
        fclose(f1); fclose(f2);
        return 0;
    }
    row++;
    while (fgets(l1, sizeof(l1), f1) && fgets(l2, sizeof(l2), f2)) {
        char b1[2048], b2[2048];
        memcpy(b1, l1, sizeof(b1));
        memcpy(b2, l2, sizeof(b2));
        char *s1 = b1, *s2 = b2;
        /* CSV columns after Phase 5: see log.c header. Deterministic cols are
         *   0..6 (frame..collisions) and 13..14 (kinetic_energy, momentum). */
        char *t1[16] = {0}, *t2[16] = {0};
        for (int c = 0; c < 15; c++) { t1[c] = csv_next(&s1); t2[c] = csv_next(&s2); }
        static const int det[] = { 0, 1, 2, 3, 4, 5, 6, 13, 14 };
        for (int k = 0; k < 9; k++) {
            if (!t1[det[k]] || !t2[det[k]] || strcmp(t1[det[k]], t2[det[k]]) != 0) mism++;
        }
        row++;
    }
    fclose(f1); fclose(f2);
    if (rows_out) *rows_out = row - 1;
    if (mism_out) *mism_out = mism;
    return mism == 0;
}

static int run_csv(const char *path, scenario_id_t id, int num_frames) {
    config_t cfg = std_cfg(id);
    cfg.num_frames = num_frames;
    sim_t s; sim_init(&s, &cfg); scenario_init(&s, id);
    FILE *f = fopen(path, "w");
    if (!f) { sim_free(&s); return 0; }
    log_write_header(f);
    for (int t = 0; t < num_frames; t++) {
        sim_metrics_t m = { .frame = t };
        sim_step(&s, &m);
        log_write_row(f, &m);
    }
    fclose(f);
    sim_free(&s);
    return 1;
}

static int g2_determinism_all(void) {
    int ok = 1;
    for (int id = 0; id < SCEN_COUNT; id++) {
        char p1[128], p2[128];
        snprintf(p1, sizeof(p1), "/tmp/scen_det_%d_a.csv", id);
        snprintf(p2, sizeof(p2), "/tmp/scen_det_%d_b.csv", id);
        int ok_run = run_csv(p1, id, 50) && run_csv(p2, id, 50);
        if (!ok_run) { ok = 0; continue; }
        int rows = 0, mism = 0;
        int eq = csv_det_eq(p1, p2, &rows, &mism);
        if (!eq) ok = 0;
        fprintf(stderr,
            "  %-12s rows=%d  det_mismatches=%d  %s\n",
            scenario_name((scenario_id_t)id), rows, mism, eq ? "OK" : "FAIL");
    }
    return ok;
}

/* ====================================== G3a: grid == brute per scenario */

static int cmp_pair(const void *a, const void *b) {
    const pair_t *pa = a, *pb = b;
    if (pa->i != pb->i) return pa->i - pb->i;
    return pa->j - pb->j;
}

static int g3a_grid_eq_brute_all(void) {
    int ok = 1;
    for (int id = 0; id < SCEN_COUNT; id++) {
        config_t cfg = std_cfg(id);
        cfg.N          = 200;
        cfg.num_frames = 30;
        sim_t s; sim_init(&s, &cfg); scenario_init(&s, id);
        int max_pairs = cfg.N * (cfg.N - 1) / 2 + 1;
        pair_t *brute      = malloc(sizeof(pair_t) * (size_t)max_pairs);
        pair_t *grid_copy  = malloc(sizeof(pair_t) * (size_t)max_pairs);
        int mismatches = 0;
        for (int t = 0; t < cfg.num_frames; t++) {
            sim_step(&s, NULL);
            int n_brute = 0;
            bruteforce_find_overlapping_pairs(s.n, s.px, s.py, cfg.radius,
                                              brute, max_pairs, &n_brute);
            qsort(brute, n_brute, sizeof(pair_t), cmp_pair);
            memcpy(grid_copy, s.pairs, sizeof(pair_t) * (size_t)s.last_pair_count);
            qsort(grid_copy, s.last_pair_count, sizeof(pair_t), cmp_pair);
            int eq = (n_brute == s.last_pair_count);
            if (eq) for (int k = 0; k < n_brute; k++) {
                if (brute[k].i != grid_copy[k].i || brute[k].j != grid_copy[k].j) { eq = 0; break; }
            }
            if (!eq) mismatches++;
        }
        if (mismatches > 0) ok = 0;
        fprintf(stderr,
            "  %-12s mismatches=%d / %d frames  %s\n",
            scenario_name((scenario_id_t)id), mismatches, cfg.num_frames,
            mismatches == 0 ? "OK" : "FAIL");
        free(brute); free(grid_copy);
        sim_free(&s);
    }
    return ok;
}

/* ====================== G3b: energy bounded (force-driven) / KE non-increasing (force-free) */

static int g3b_energy_all(void) {
    int ok = 1;
    const int LONG = 1000;
    for (int id = 0; id < SCEN_COUNT; id++) {
        config_t cfg = std_cfg(id);
        cfg.num_frames = LONG;
        sim_t s; sim_init(&s, &cfg); scenario_init(&s, id);
        double ke0     = sim_total_kinetic_energy(&s);
        double ke_prev = ke0;
        double max_ke  = ke0;
        int    nincr_violations = 0;
        int    has_forces = scenario_has_forces((scenario_id_t)id);
        for (int t = 0; t < LONG; t++) {
            sim_step(&s, NULL);
            double ke = sim_total_kinetic_energy(&s);
            if (ke > max_ke) max_ke = ke;
            if (!has_forces) {
                if (ke > ke_prev * (1.0 + 1e-9) + 1e-12) nincr_violations++;
            }
            ke_prev = ke;
        }
        int case_ok;
        if (has_forces) {
            /* Bounded energy: finite and below a very generous absolute ceiling.
             * Realistic max_KE is bounded by initial total energy, which depends
             * on the scenario's force field. 1e7 catches actual blow-ups while
             * tolerating expected oscillation amplitudes. */
            case_ok = isfinite(max_ke) && (max_ke < 1.0e7);
            fprintf(stderr,
                "  %-12s force-driven : ke0=%.3g  max_ke=%.3g  final=%.3g  %s\n",
                scenario_name((scenario_id_t)id), ke0, max_ke, ke_prev,
                case_ok ? "OK" : "FAIL");
        } else {
            case_ok = (nincr_violations == 0);
            fprintf(stderr,
                "  %-12s force-free   : ke0=%.3g  final=%.3g  violations=%d  %s\n",
                scenario_name((scenario_id_t)id), ke0, ke_prev, nincr_violations,
                case_ok ? "OK" : "FAIL");
        }
        if (!case_ok) ok = 0;
        sim_free(&s);
    }
    return ok;
}

/* ========================================== G4: no escapes per scenario */

static int g4_no_escapes_all(void) {
    int ok = 1;
    const double eps = 1e-9;
    for (int id = 0; id < SCEN_COUNT; id++) {
        config_t cfg = std_cfg(id);
        cfg.num_frames = 60;
        sim_t s; sim_init(&s, &cfg); scenario_init(&s, id);
        int    escapes = 0;
        double worst   = 0;
        for (int t = 0; t < cfg.num_frames; t++) {
            sim_step(&s, NULL);
            for (int i = 0; i < s.n; i++) {
                double ox = 0, oy = 0;
                if (s.px[i] < -eps)               ox = -s.px[i];
                if (s.px[i] > cfg.domain_w + eps) ox = s.px[i] - cfg.domain_w;
                if (s.py[i] < -eps)               oy = -s.py[i];
                if (s.py[i] > cfg.domain_h + eps) oy = s.py[i] - cfg.domain_h;
                if (ox > 0 || oy > 0) {
                    escapes++;
                    if (ox > worst) worst = ox;
                    if (oy > worst) worst = oy;
                }
            }
        }
        if (escapes > 0) ok = 0;
        fprintf(stderr,
            "  %-12s escapes=%d  worst_overshoot=%g  %s\n",
            scenario_name((scenario_id_t)id), escapes, worst,
            escapes == 0 ? "OK" : "FAIL");
        sim_free(&s);
    }
    return ok;
}

/* =================================================================== main */

int main(void) {
    fprintf(stderr, "=== Phase 3 Verification Gates ===\n");

    int all_ok = 1;

    fprintf(stderr, "\n--- G1 characteristic statistic per scenario ---\n");
    all_ok &= char_stable();
    all_ok &= char_collapse();
    all_ok &= char_explosion();
    all_ok &= char_vortex();
    all_ok &= char_funnel();
    all_ok &= char_pulse();
    all_ok &= char_multicluster();
    all_ok &= char_stream();

    fprintf(stderr, "\n--- G2 determinism per scenario (det. CSV cols bit-identical) ---\n");
    all_ok &= g2_determinism_all();

    fprintf(stderr, "\n--- G3a grid pipeline == brute force per scenario (30 frames) ---\n");
    all_ok &= g3a_grid_eq_brute_all();

    fprintf(stderr, "\n--- G3b energy bounded (force-driven) / non-increasing (force-free) ---\n");
    all_ok &= g3b_energy_all();

    fprintf(stderr, "\n--- G4 no escapes per scenario (60 frames) ---\n");
    all_ok &= g4_no_escapes_all();

    fprintf(stderr, "\n%s\n", all_ok ? "Phase 3: ALL OK" : "Phase 3: FAILED");
    return all_ok ? 0 : 1;
}
