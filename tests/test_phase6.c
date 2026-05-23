/*
 * Phase 6 verification gates (spec §10 build phase 6 -- the 3D extension).
 *
 *   G1. The 3D grid pipeline (counting-sort + 13-neighbor half-stencil broad
 *       phase + sphere narrow phase) returns EXACTLY the brute-force
 *       overlapping-pair set on small N -- no missed, no spurious pairs --
 *       across several cell sizes (ell >= 2r) and several frames.
 *
 *   G2. 3D collision physics: an isolated pair conserves total momentum, and
 *       at e=1 conserves kinetic energy; over a force-free run global KE is
 *       non-increasing for e <= 1.
 *
 *   G3. The closed-form C output for d=3 matches math_kernel_reference.py's
 *       d=3 canonical vectors to <= 1e-12 relative.
 *
 *   G4. The 3D method runs end-to-end on all 3D scenarios (stable, collapse,
 *       explosion, pulse, settling), deterministically, with ell bounded in
 *       [ell_min, ell_max] and no runaway, in Modes A and B.
 *
 *   G5. 3D cost-model calibration: fit T = c0 + a*M + b*S over an ell sweep
 *       across representative 3D frames; report R^2 and the fitted a, b, c0.
 */

#define _POSIX_C_SOURCE 199309L
#include "adapter.h"
#include "config.h"
#include "sim.h"
#include "scenarios.h"
#include "grid.h"
#include "oracle.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static config_t cfg3d(scenario_id_t id, int N, double box, int frames) {
    config_t cfg; config_default(&cfg);
    cfg.dim         = 3;
    cfg.N           = N;
    cfg.domain_w    = box;
    cfg.domain_h    = box;
    cfg.domain_d    = box;
    cfg.radius      = 0.5;
    cfg.cell_size   = 2.0;
    cfg.dt          = 1.0 / 60.0;
    cfg.init_speed  = 10.0;
    cfg.restitution = 0.95;
    cfg.num_frames  = frames;
    cfg.seed        = 31;
    snprintf(cfg.scenario, sizeof(cfg.scenario), "%s", scenario_name(id));
    return cfg;
}

/* ===================================================== G1 grid == brute */

static int cmp_pair(const void *a, const void *b) {
    const pair_t *pa = (const pair_t *)a, *pb = (const pair_t *)b;
    if (pa->i != pb->i) return pa->i - pb->i;
    return pa->j - pb->j;
}

static int pairsets_equal(pair_t *A, int na, pair_t *B, int nb) {
    if (na != nb) return 0;
    qsort(A, na, sizeof(pair_t), cmp_pair);
    qsort(B, nb, sizeof(pair_t), cmp_pair);
    for (int k = 0; k < na; k++)
        if (A[k].i != B[k].i || A[k].j != B[k].j) return 0;
    return 1;
}

static int g1_grid_vs_brute(void) {
    /* Dense-ish small 3D box so collisions create genuine overlaps. */
    config_t cfg = cfg3d(SCEN_COLLAPSE, 1500, 40.0, 40);
    sim_t s; sim_init(&s, &cfg);
    scenario_init(&s, SCEN_COLLAPSE);

    const int cap = 64 * cfg.N;
    pair_t *grid_pairs  = (pair_t *)malloc(sizeof(pair_t) * (size_t)cap);
    pair_t *brute_pairs = (pair_t *)malloc(sizeof(pair_t) * (size_t)cap);
    double ells[3] = { 1.0, 2.0, 4.0 };     /* all >= 2r = 1.0 */

    int ok = 1, total_overlaps = 0, comparisons = 0, mismatches = 0;
    for (int t = 0; t < cfg.num_frames; t++) {
        sim_step(&s, NULL);
        if (t % 8 != 0) continue;            /* sample a handful of frames */
        for (int e = 0; e < 3; e++) {
            grid_t g;
            grid_init3d(&g, s.n, cfg.domain_w, cfg.domain_h, cfg.domain_d, ells[e]);
            grid_build3d(&g, s.px, s.py, s.pz);
            int ng = 0;
            if (grid_find_overlapping_pairs3d(&g, s.px, s.py, s.pz, cfg.radius,
                                              grid_pairs, cap, &ng) != 0) {
                fprintf(stderr, "  G1: grid pair buffer overflow\n"); ok = 0;
            }
            int nb = 0;
            bruteforce_find_overlapping_pairs3d(s.n, s.px, s.py, s.pz, cfg.radius,
                                                brute_pairs, cap, &nb);
            if (!pairsets_equal(grid_pairs, ng, brute_pairs, nb)) {
                fprintf(stderr, "  G1 MISMATCH frame=%d ell=%.1f grid=%d brute=%d\n",
                        t, ells[e], ng, nb);
                mismatches++; ok = 0;
            }
            total_overlaps += nb;
            comparisons++;
            grid_free(&g);
        }
    }
    /* Guard against a vacuous pass (the scenario must actually produce overlaps). */
    if (total_overlaps == 0) { fprintf(stderr, "  G1: no overlaps seen (vacuous)\n"); ok = 0; }

    fprintf(stderr,
        "  G1 grid==brute (3D): %d comparisons over ells{1,2,4}, %d overlaps total, "
        "%d mismatches  %s\n",
        comparisons, total_overlaps, mismatches, ok ? "OK" : "FAIL");
    free(grid_pairs); free(brute_pairs);
    sim_free(&s);
    return ok;
}

/* ===================================================== G2 physics */

static int g2_physics(void) {
    int ok = 1;

    /* (a) isolated 3D pair: momentum conserved; e=1 => KE conserved. */
    double worst_p = 0, worst_ke = 0;
    for (int trial = 0; trial < 200; trial++) {
        rng_state_t rng; rng_seed(&rng, 1000 + trial);
        double px[2], py[2], pz[2], vx[2], vy[2], vz[2];
        /* two particles within 2r, approaching */
        px[0]=0; py[0]=0; pz[0]=0;
        double th = rng_uniform(&rng,0,6.283), u = rng_uniform(&rng,-1,1);
        double sp = sqrt(1-u*u), d = 0.8;     /* < 2r = 1 -> overlapping */
        px[1]=d*sp*cos(th); py[1]=d*sp*sin(th); pz[1]=d*u;
        for (int q=0;q<2;q++){ vx[q]=rng_uniform(&rng,-5,5); vy[q]=rng_uniform(&rng,-5,5); vz[q]=rng_uniform(&rng,-5,5);}
        double px_tot0 = vx[0]+vx[1], py_tot0 = vy[0]+vy[1], pz_tot0 = vz[0]+vz[1];
        double ke0 = 0.5*(vx[0]*vx[0]+vy[0]*vy[0]+vz[0]*vz[0]
                         +vx[1]*vx[1]+vy[1]*vy[1]+vz[1]*vz[1]);
        pair_t pr = { 0, 1 };
        resolve_pair_collisions3d(vx,vy,vz, px,py,pz, &pr, 1, 1.0);   /* e=1 */
        double dpx = (vx[0]+vx[1])-px_tot0, dpy=(vy[0]+vy[1])-py_tot0, dpz=(vz[0]+vz[1])-pz_tot0;
        double dp = sqrt(dpx*dpx+dpy*dpy+dpz*dpz);
        double ke1 = 0.5*(vx[0]*vx[0]+vy[0]*vy[0]+vz[0]*vz[0]
                         +vx[1]*vx[1]+vy[1]*vy[1]+vz[1]*vz[1]);
        double rel_ke = fabs(ke1-ke0)/(ke0+1e-30);
        if (dp > worst_p) worst_p = dp;
        if (rel_ke > worst_ke) worst_ke = rel_ke;
    }
    int pair_ok = (worst_p < 1e-12) && (worst_ke < 1e-12);
    fprintf(stderr, "  G2a isolated pair: worst |dP|=%.2e (momentum), worst dKE/KE=%.2e (e=1)  %s\n",
            worst_p, worst_ke, pair_ok ? "OK" : "FAIL");
    if (!pair_ok) ok = 0;

    /* (b) force-free 3D run, e<=1 => global KE non-increasing. */
    config_t cfg = cfg3d(SCEN_STABLE, 3000, 60.0, 150);
    cfg.restitution = 0.95;
    sim_t s; sim_init(&s, &cfg);
    scenario_init(&s, SCEN_STABLE);
    double ke_prev = sim_total_kinetic_energy(&s), ke0 = ke_prev;
    int violations = 0;
    for (int t = 0; t < cfg.num_frames; t++) {
        sim_step(&s, NULL);
        double ke = sim_total_kinetic_energy(&s);
        if (ke > ke_prev * (1.0 + 1e-9) + 1e-12) violations++;
        ke_prev = ke;
    }
    int run_ok = (violations == 0);
    fprintf(stderr, "  G2b force-free run (e=0.95): ke0=%.4g final=%.4g  KE-increase frames=%d  %s\n",
            ke0, ke_prev, violations, run_ok ? "OK" : "FAIL");
    if (!run_ok) ok = 0;
    sim_free(&s);
    return ok;
}

/* ===================================================== G3 d=3 kernel */

static int g3_kernel_d3(void) {
    /* d=3 canonical vectors from math_kernel_reference.py (full precision). */
    struct { int d; double D2, ell, B, Q, expected; } cases[] = {
        { 3, 3.0, 1.0,  9.0,  1.0, 1.4422495703074083 },
        { 3, 1.0, 4.0,  8.0, 50.0, 3.329433160230254  },
        { 3, 2.0, 2.0,  5.0,  5.0, 2.168943542395397  },
    };
    int n = (int)(sizeof(cases)/sizeof(cases[0])), ok = 1;
    double worst = 0;
    for (int k = 0; k < n; k++) {
        double got = adapter_optimal_ell(cases[k].ell, cases[k].d, cases[k].D2,
                                         cases[k].B, cases[k].Q);
        double rel = fabs(got - cases[k].expected) / cases[k].expected;
        if (rel > worst) worst = rel;
        if (rel > 1e-12) ok = 0;
    }
    fprintf(stderr, "  G3 d=3 kernel vs reference: worst rel err = %.2e  %s\n",
            worst, ok ? "OK" : "FAIL");
    return ok;
}

/* ===================================================== G4 end-to-end */

static int run3d_adapter(scenario_id_t id, adapter_mode_t mode,
                         double *ell_lo, double *ell_hi, int *oor, double *mean_oh_pct) {
    config_t cfg = cfg3d(id, 3000, 100.0, 80);
    sim_t s; sim_init(&s, &cfg);
    scenario_init(&s, id);
    adapter_t a;
    adapter_init_default(&a, mode, 3, cfg.domain_w, cfg.domain_h, cfg.cell_size, cfg.radius);
    a.domain_d = cfg.domain_d;

    double lo = a.ell_cur, hi = a.ell_cur, sum_oh = 0, sum_td = 0;
    int bad = 0;
    for (int t = 0; t < cfg.num_frames; t++) {
        sim_metrics_t m = { .frame = t };
        sim_step(&s, &m);
        double t_M = m.t_M + m.t_broad_iter;
        double t_S = (m.t_broad - m.t_broad_iter) + m.t_narrow;
        if (t_S < 0) t_S = 0;
        sum_td += (m.t_grid + m.t_broad + m.t_narrow);
        double ell_next = adapter_step(&a, t_M, t_S, m.S, m.num_cells,
                                       s.n, s.px, s.py, s.pz);
        sum_oh += a.last_overhead_s;
        if (!isfinite(ell_next) || ell_next < a.ell_min - 1e-9 || ell_next > a.ell_max + 1e-9) bad++;
        if (a.ell_cur < lo) lo = a.ell_cur;
        if (a.ell_cur > hi) hi = a.ell_cur;
        if (a.last_resized) grid_resize(&s.grid, ell_next);
    }
    *ell_lo = lo; *ell_hi = hi; *oor = bad;
    *mean_oh_pct = 100.0 * (sum_oh / cfg.num_frames) / (sum_td / cfg.num_frames);
    sim_free(&s);
    return bad;
}

static int g4_end_to_end(void) {
    scenario_id_t ids[] = { SCEN_STABLE, SCEN_COLLAPSE, SCEN_EXPLOSION, SCEN_PULSE, SCEN_SETTLING };
    int ok = 1;
    for (int k = 0; k < (int)(sizeof(ids)/sizeof(ids[0])); k++) {
        double la, ha, lb, hb, oha, ohb; int oora, oorb;
        run3d_adapter(ids[k], ADAPTER_MODE_A, &la, &ha, &oora, &oha);
        run3d_adapter(ids[k], ADAPTER_MODE_B, &lb, &hb, &oorb, &ohb);
        int case_ok = (oora == 0) && (oorb == 0)
                   && (ha/la < 1e3) && (hb/lb < 1e3);
        fprintf(stderr,
            "  %-10s : modeA ell[%.2f,%.2f] oh=%.2f%%  |  modeB ell[%.2f,%.2f] oh=%.2f%%  %s\n",
            scenario_name(ids[k]), la, ha, oha, lb, hb, ohb, case_ok ? "OK" : "FAIL");
        if (!case_ok) ok = 0;
    }
    return ok;
}

/* ===================================================== G5 3D calibration */

static int g5_calibration(void) {
    /* Fit T = c0 + a*M + b*S over an ell sweep across a few 3D frames. */
    config_t cfg = cfg3d(SCEN_COLLAPSE, 4000, 100.0, 60);
    double ells[32];
    int n_ells = oracle_make_geo_sweep(2.0*cfg.radius, 25.0, 1.3, ells, 32);

    sim_t s; sim_init(&s, &cfg);
    scenario_init(&s, SCEN_COLLAPSE);
    int n_frames = 6, K = 5;
    int cap = n_frames * n_ells, ns = 0;
    double *M = malloc(sizeof(double)*cap), *S = malloc(sizeof(double)*cap), *T = malloc(sizeof(double)*cap);
    oracle_point_t *pts = malloc(sizeof(*pts) * n_ells);
    int targets[6];
    for (int k=0;k<n_frames;k++) targets[k] = (int)((double)(k+1)*cfg.num_frames/(n_frames+1)+0.5);
    int tk=0;
    for (int t=0;t<cfg.num_frames && tk<n_frames;t++){
        sim_step(&s, NULL);
        if (t != targets[tk]) continue;
        oracle_eval_sweep3d(s.n, s.px, s.py, s.pz, cfg.radius,
                            cfg.domain_w, cfg.domain_h, cfg.domain_d, ells, n_ells, K, pts);
        for (int i=0;i<n_ells && ns<cap;i++){ M[ns]=pts[i].num_cells; S[ns]=(double)pts[i].S; T[ns]=pts[i].t_detect; ns++; }
        tk++;
    }
    cost_fit_t f = adapter_fit_cost_model(M, S, T, ns);
    int ok = f.ok && f.a > 0.0 && f.b > 0.0 && f.r2 > 0.90;
    fprintf(stderr,
        "  G5 3D calibration: a=%.3e b=%.3e c0=%.3e  R^2=%.4f  RMSE=%.2fus  mean|r|=%.1f%%  [n=%d]  %s\n",
        f.a, f.b, f.c0, f.r2, f.rmse*1e6, f.mean_rel_resid*100.0, f.n, ok ? "OK" : "FAIL");
    free(M); free(S); free(T); free(pts);
    sim_free(&s);
    return ok;
}

/* =================================================================== main */

int main(void) {
    fprintf(stderr, "=== Phase 6 Verification Gates (3D) ===\n");
    int all_ok = 1;

    fprintf(stderr, "\n--- G1 3D grid pipeline == brute force ---\n");
    all_ok &= g1_grid_vs_brute();

    fprintf(stderr, "\n--- G2 3D collision physics ---\n");
    all_ok &= g2_physics();

    fprintf(stderr, "\n--- G3 closed-form d=3 vs reference ---\n");
    all_ok &= g3_kernel_d3();

    fprintf(stderr, "\n--- G4 3D end-to-end on all 3D scenarios ---\n");
    all_ok &= g4_end_to_end();

    fprintf(stderr, "\n--- G5 3D cost-model calibration ---\n");
    all_ok &= g5_calibration();

    fprintf(stderr, "\n%s\n", all_ok ? "Phase 6: ALL OK" : "Phase 6: FAILED");
    return all_ok ? 0 : 1;
}
