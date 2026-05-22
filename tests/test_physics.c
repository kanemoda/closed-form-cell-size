/*
 * Verification gate 2 (spec section 10):
 *   "Collision physics is correct: a pairwise collision conserves total
 *    momentum and (at restitution=1) total kinetic energy; global kinetic
 *    energy is non-increasing for restitution <= 1."
 *
 * Tests:
 *   T1.  head-on pair, e=1   : KE & P conserved; velocities swap
 *   T2.  head-on pair, e=0.5 : P conserved; KE_after/KE_before = 0.25
 *   T3.  oblique  pair, e=1  : KE & P conserved
 *   T4.  separating pair     : no impulse applied
 *   T5.  isolated pair via full sim_step (no walls hit): KE & P conserved
 *   T6.  multi-particle interior (no walls hit), e=1: KE & P conserved
 *   T7.  full sim with walls, e=0.5: KE non-increasing per frame
 */

#include "rng.h"
#include "config.h"
#include "sim.h"
#include "grid.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int approx_eq_abs(double a, double b, double tol) {
    return fabs(a - b) <= tol;
}

/* ------------------------------- T1: head-on pair, e=1 */
static int t1_pair_e1_headon(void) {
    double px[2] = { 10.0, 10.99 };   /* dist 0.99 < 2r=1 -> overlapping */
    double py[2] = { 10.0, 10.0  };
    double vx[2] = {  1.0, -1.0  };
    double vy[2] = {  0.0,  0.0  };

    double ke0 = 0.5 * (vx[0]*vx[0] + vy[0]*vy[0] + vx[1]*vx[1] + vy[1]*vy[1]);
    double Px0 = vx[0] + vx[1], Py0 = vy[0] + vy[1];

    pair_t p = { 0, 1 };
    resolve_pair_collisions(vx, vy, px, py, &p, 1, 1.0);

    double ke1 = 0.5 * (vx[0]*vx[0] + vy[0]*vy[0] + vx[1]*vx[1] + vy[1]*vy[1]);
    double Px1 = vx[0] + vx[1], Py1 = vy[0] + vy[1];

    int ok = 1;
    ok &= approx_eq_abs(ke1, ke0, 1e-14);
    ok &= approx_eq_abs(Px1, Px0, 1e-14);
    ok &= approx_eq_abs(Py1, Py0, 1e-14);
    /* head-on, equal mass, e=1: velocities swap */
    ok &= approx_eq_abs(vx[0], -1.0, 1e-14);
    ok &= approx_eq_abs(vx[1],  1.0, 1e-14);

    fprintf(stderr, "  T1 head-on e=1: KE %.6g->%.6g  P (%g,%g)->(%g,%g)  v0->%g v1->%g  %s\n",
            ke0, ke1, Px0, Py0, Px1, Py1, vx[0], vx[1], ok ? "OK" : "FAIL");
    return ok;
}

/* --------------------------- T2: head-on pair, e=0.5 -- KE ratio = e^2 = 0.25 */
static int t2_pair_e_inelastic(void) {
    double px[2] = { 10.0, 10.99 };
    double py[2] = { 10.0, 10.0  };
    double vx[2] = {  1.0, -1.0  };
    double vy[2] = {  0.0,  0.0  };

    double ke0 = 0.5 * (vx[0]*vx[0] + vy[0]*vy[0] + vx[1]*vx[1] + vy[1]*vy[1]);
    double Px0 = vx[0] + vx[1], Py0 = vy[0] + vy[1];

    pair_t p = { 0, 1 };
    resolve_pair_collisions(vx, vy, px, py, &p, 1, 0.5);

    double ke1 = 0.5 * (vx[0]*vx[0] + vy[0]*vy[0] + vx[1]*vx[1] + vy[1]*vy[1]);
    double Px1 = vx[0] + vx[1], Py1 = vy[0] + vy[1];

    int ok = 1;
    ok &= approx_eq_abs(Px1, Px0, 1e-14);
    ok &= approx_eq_abs(Py1, Py0, 1e-14);
    /* For zero-momentum head-on: KE_after/KE_before = e^2 = 0.25 */
    ok &= approx_eq_abs(ke1 / ke0, 0.25, 1e-14);

    fprintf(stderr, "  T2 head-on e=0.5: KE %.6g->%.6g (ratio=%.6f)  P (%g,%g)->(%g,%g)  %s\n",
            ke0, ke1, ke1 / ke0, Px0, Py0, Px1, Py1, ok ? "OK" : "FAIL");
    return ok;
}

/* ------------------------- T3: oblique pair, e=1 */
static int t3_pair_e1_oblique(void) {
    double px[2] = { 0.0, 0.7 };       /* dist sqrt(0.98) ~ 0.99 -> overlapping */
    double py[2] = { 0.0, 0.7 };
    double vx[2] = { 1.0, 0.0 };
    double vy[2] = { 0.5, 0.0 };

    double ke0 = 0.5 * (vx[0]*vx[0] + vy[0]*vy[0] + vx[1]*vx[1] + vy[1]*vy[1]);
    double Px0 = vx[0] + vx[1], Py0 = vy[0] + vy[1];

    pair_t p = { 0, 1 };
    resolve_pair_collisions(vx, vy, px, py, &p, 1, 1.0);

    double ke1 = 0.5 * (vx[0]*vx[0] + vy[0]*vy[0] + vx[1]*vx[1] + vy[1]*vy[1]);
    double Px1 = vx[0] + vx[1], Py1 = vy[0] + vy[1];

    int ok = 1;
    ok &= approx_eq_abs(ke1, ke0, 1e-14);
    ok &= approx_eq_abs(Px1, Px0, 1e-14);
    ok &= approx_eq_abs(Py1, Py0, 1e-14);

    fprintf(stderr, "  T3 oblique e=1: KE %.6g->%.6g  P (%g,%g)->(%g,%g)  %s\n",
            ke0, ke1, Px0, Py0, Px1, Py1, ok ? "OK" : "FAIL");
    return ok;
}

/* ----------------- T4: pair already separating -- impulse skipped */
static int t4_pair_separating(void) {
    double px[2] = { 10.0, 10.99 };
    double py[2] = { 10.0, 10.0  };
    double vx[2] = { -1.0,  1.0  };   /* moving apart along normal */
    double vy[2] = {  0.0,  0.0  };

    double ke0 = 0.5 * (vx[0]*vx[0] + vy[0]*vy[0] + vx[1]*vx[1] + vy[1]*vy[1]);

    pair_t p = { 0, 1 };
    resolve_pair_collisions(vx, vy, px, py, &p, 1, 1.0);

    double ke1 = 0.5 * (vx[0]*vx[0] + vy[0]*vy[0] + vx[1]*vx[1] + vy[1]*vy[1]);
    int ok = approx_eq_abs(vx[0], -1.0, 1e-14) && approx_eq_abs(vx[1], 1.0, 1e-14)
          && approx_eq_abs(ke1, ke0, 1e-14);
    fprintf(stderr, "  T4 separating: no impulse  %s\n", ok ? "OK" : "FAIL");
    return ok;
}

/* ------- T5: isolated pair through sim_step, no walls hit */
static int t5_isolated_via_simstep(void) {
    config_t cfg;
    config_default(&cfg);
    cfg.N           = 2;
    cfg.domain_w    = 200;
    cfg.domain_h    = 200;
    cfg.radius      = 0.5;
    cfg.cell_size   = 2.0;
    cfg.dt          = 0.01;
    cfg.init_speed  = 1.0;
    cfg.num_frames  = 80;     /* finishes well inside the domain */
    cfg.seed        = 1;
    cfg.restitution = 1.0;

    sim_t s;
    sim_init(&s, &cfg);
    /* place manually -- bypasses sim_init_random */
    s.px[0] = 99.5;   s.py[0] = 100.0;  s.vx[0] =  1.0;  s.vy[0] = 0.0;
    s.px[1] = 100.5;  s.py[1] = 100.0;  s.vx[1] = -1.0;  s.vy[1] = 0.0;

    double ke0 = sim_total_kinetic_energy(&s);
    double Px0, Py0; sim_total_momentum(&s, &Px0, &Py0);

    int collided = 0;
    for (int t = 0; t < cfg.num_frames; t++) {
        sim_step(&s, NULL);
        if (s.last_pair_count > 0) collided = 1;
    }

    double ke1 = sim_total_kinetic_energy(&s);
    double Px1, Py1; sim_total_momentum(&s, &Px1, &Py1);

    int ok = 1;
    ok &= collided;                                    /* sanity: collision must have happened */
    ok &= approx_eq_abs(ke1, ke0,        1e-12);
    ok &= approx_eq_abs(Px1, Px0,        1e-12);
    ok &= approx_eq_abs(Py1, Py0,        1e-12);
    /* and verify no walls were hit: particles in interior */
    for (int i = 0; i < s.n; i++) {
        if (s.px[i] < 2.0 * cfg.radius || s.px[i] > cfg.domain_w - 2.0 * cfg.radius ||
            s.py[i] < 2.0 * cfg.radius || s.py[i] > cfg.domain_h - 2.0 * cfg.radius) {
            ok = 0;
        }
    }
    fprintf(stderr, "  T5 isolated-pair via sim_step e=1: collided=%d KE %.6g->%.6g  P (%g,%g)->(%g,%g)  %s\n",
            collided, ke0, ke1, Px0, Py0, Px1, Py1, ok ? "OK" : "FAIL");

    sim_free(&s);
    return ok;
}

/* ------- T6: many particles in interior, no walls hit, e=1 -- KE + P conserved */
static int t6_global_no_walls_e1(void) {
    /* Constrain particles to [200, 280]x[200, 280], domain 500x500, init_speed=1,
     * dt=0.01, num_frames=50: max travel is 0.5 -> nowhere near walls at x=0,500. */
    config_t cfg;
    config_default(&cfg);
    cfg.N           = 200;
    cfg.domain_w    = 500;
    cfg.domain_h    = 500;
    cfg.radius      = 0.5;
    cfg.cell_size   = 2.0;
    cfg.dt          = 0.01;
    cfg.init_speed  = 1.0;
    cfg.num_frames  = 50;
    cfg.seed        = 77;
    cfg.restitution = 1.0;

    sim_t s;
    sim_init(&s, &cfg);

    /* Manual placement in central window, with rejection sampling. */
    rng_state_t rng;
    rng_seed(&rng, 77);
    double rr = (2 * cfg.radius) * (2 * cfg.radius);
    int placed = 0;
    while (placed < cfg.N) {
        double x = rng_uniform(&rng, 200.0, 280.0);
        double y = rng_uniform(&rng, 200.0, 280.0);
        int ok = 1;
        for (int j = 0; j < placed; j++) {
            double dx = x - s.px[j], dy = y - s.py[j];
            if (dx * dx + dy * dy < rr) { ok = 0; break; }
        }
        if (ok) {
            s.px[placed] = x;
            s.py[placed] = y;
            double th = rng_uniform(&rng, 0.0, 2.0 * M_PI);
            s.vx[placed] = cfg.init_speed * cos(th);
            s.vy[placed] = cfg.init_speed * sin(th);
            placed++;
        }
    }

    double ke0 = sim_total_kinetic_energy(&s);
    double Px0, Py0; sim_total_momentum(&s, &Px0, &Py0);

    long total_pairs = 0;
    for (int t = 0; t < cfg.num_frames; t++) {
        sim_step(&s, NULL);
        total_pairs += s.last_pair_count;
    }

    double ke1 = sim_total_kinetic_energy(&s);
    double Px1, Py1; sim_total_momentum(&s, &Px1, &Py1);

    int hit_wall = 0;
    for (int i = 0; i < s.n; i++) {
        if (s.px[i] < 2.0 || s.px[i] > cfg.domain_w - 2.0 ||
            s.py[i] < 2.0 || s.py[i] > cfg.domain_h - 2.0) {
            hit_wall = 1;
        }
    }

    int ok = 1;
    ok &= !hit_wall;
    /* Tolerance scaled by KE and #collisions: each collision adds ~eps*KE
     * worth of round-off; accumulating linearly. With ~total_pairs=hundreds,
     * 1e-10 absolute is generous. */
    double ke_drift = fabs(ke1 - ke0);
    double p_drift  = fabs(Px1 - Px0) + fabs(Py1 - Py0);
    ok &= ke_drift < 1e-9 * (ke0 + 1.0);
    ok &= p_drift  < 1e-9 * (fabs(Px0) + fabs(Py0) + 1.0);

    fprintf(stderr,
        "  T6 N=%d interior e=1 (no walls hit=%d): total_pairs=%ld KE drift=%.3g P drift=%.3g  %s\n",
        cfg.N, hit_wall, total_pairs, ke_drift, p_drift, ok ? "OK" : "FAIL");

    sim_free(&s);
    return ok;
}

/* ------- T7: full sim with walls, e=0.5 -- KE strictly non-increasing per frame */
static int t7_global_ke_nonincreasing(void) {
    config_t cfg;
    config_default(&cfg);
    cfg.N           = 500;
    cfg.domain_w    = 80;
    cfg.domain_h    = 80;
    cfg.radius      = 0.5;
    cfg.cell_size   = 2.0;
    cfg.dt          = 0.01;
    cfg.init_speed  = 10.0;
    cfg.num_frames  = 300;
    cfg.seed        = 99;
    cfg.restitution = 0.5;

    sim_t s;
    sim_init(&s, &cfg);
    sim_init_random(&s);

    double ke0       = sim_total_kinetic_energy(&s);
    double ke_prev   = ke0;
    int    violations = 0;
    double worst_up  = 0.0;
    long   total_pairs = 0;
    for (int t = 0; t < cfg.num_frames; t++) {
        sim_step(&s, NULL);
        total_pairs += s.last_pair_count;
        double ke   = sim_total_kinetic_energy(&s);
        double diff = ke - ke_prev;
        /* allow noise on the order of round-off times KE */
        if (diff > 1e-12 * (ke0 + 1.0)) {
            violations++;
            if (diff > worst_up) worst_up = diff;
        }
        ke_prev = ke;
    }
    double ke_final = sim_total_kinetic_energy(&s);

    int ok = (violations == 0) && (ke_final <= ke0 + 1e-12 * ke0);
    fprintf(stderr,
        "  T7 N=%d e=0.5 walls %d frames: total_pairs=%ld KE %.6g -> %.6g (ratio=%.4f) violations=%d worst_up=%.3g  %s\n",
        cfg.N, cfg.num_frames, total_pairs,
        ke0, ke_final, ke_final / ke0, violations, worst_up,
        ok ? "OK" : "FAIL");

    sim_free(&s);
    return ok;
}

int main(void) {
    fprintf(stderr, "=== Gate 2: collision physics ===\n");

    int all_ok = 1;
    all_ok &= t1_pair_e1_headon();
    all_ok &= t2_pair_e_inelastic();
    all_ok &= t3_pair_e1_oblique();
    all_ok &= t4_pair_separating();
    all_ok &= t5_isolated_via_simstep();
    all_ok &= t6_global_no_walls_e1();
    all_ok &= t7_global_ke_nonincreasing();

    fprintf(stderr, "\n%s\n", all_ok ? "Gate 2: ALL OK" : "Gate 2: FAILED");
    return all_ok ? 0 : 1;
}
