#define _POSIX_C_SOURCE 199309L
#include "sim.h"
#include "scenarios.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static inline double ts_diff_s(const struct timespec *a, const struct timespec *b) {
    return (double)(b->tv_sec - a->tv_sec) + (double)(b->tv_nsec - a->tv_nsec) * 1e-9;
}

void sim_init(sim_t *s, const config_t *cfg) {
    s->cfg = *cfg;
    s->n   = cfg->N;
    s->px  = (double *)calloc((size_t)s->n, sizeof(double));
    s->py  = (double *)calloc((size_t)s->n, sizeof(double));
    s->vx  = (double *)calloc((size_t)s->n, sizeof(double));
    s->vy  = (double *)calloc((size_t)s->n, sizeof(double));
    if (cfg->dim == 3) {
        s->pz = (double *)calloc((size_t)s->n, sizeof(double));
        s->vz = (double *)calloc((size_t)s->n, sizeof(double));
    } else {
        s->pz = NULL;
        s->vz = NULL;
    }

    if (cfg->cell_size < 2.0 * cfg->radius) {
        fprintf(stderr,
                "sim_init: cell_size=%g < 2*radius=%g -- 1-ring stencil cannot guarantee completeness.\n",
                cfg->cell_size, 2.0 * cfg->radius);
    }

    rng_seed(&s->rng, cfg->seed);
    if (cfg->dim == 3)
        grid_init3d(&s->grid, s->n, cfg->domain_w, cfg->domain_h,
                    cfg->domain_d, cfg->cell_size);
    else
        grid_init(&s->grid, s->n, cfg->domain_w, cfg->domain_h, cfg->cell_size);

    s->pairs_cap        = 16 * s->n;
    if (s->pairs_cap < 1024) s->pairs_cap = 1024;
    s->pairs            = (pair_t *)malloc(sizeof(pair_t) * (size_t)s->pairs_cap);
    s->last_cand_count  = 0;
    s->last_pair_count  = 0;

    /* Default to the force-free 'stable' scenario; scenario_init() may
     * override before the first sim_step(). */
    s->scenario_id      = SCEN_STABLE;
    s->frame_index      = 0;
}

void sim_free(sim_t *s) {
    free(s->px);
    free(s->py);
    free(s->vx);
    free(s->vy);
    free(s->pz);
    free(s->vz);
    free(s->pairs);
    grid_free(&s->grid);
    s->px = s->py = s->vx = s->vy = s->pz = s->vz = NULL;
    s->pairs = NULL;
}

/* Legacy entry point used by Phase 1/2 tests: identical to the 'stable'
 * scenario (uniform random placement, random velocity directions). */
void sim_init_random(sim_t *s) {
    scenario_init(s, SCEN_STABLE);
}

static void wall_reflect(sim_t *s) {
    const double r = s->cfg.radius;
    const double W = s->cfg.domain_w;
    const double H = s->cfg.domain_h;
    for (int i = 0; i < s->n; i++) {
        if (s->px[i] < r)        { s->px[i] = 2.0 * r          - s->px[i]; s->vx[i] = -s->vx[i]; }
        if (s->px[i] > W - r)    { s->px[i] = 2.0 * (W - r)    - s->px[i]; s->vx[i] = -s->vx[i]; }
        if (s->py[i] < r)        { s->py[i] = 2.0 * r          - s->py[i]; s->vy[i] = -s->vy[i]; }
        if (s->py[i] > H - r)    { s->py[i] = 2.0 * (H - r)    - s->py[i]; s->vy[i] = -s->vy[i]; }
    }
}

static void wall_reflect3d(sim_t *s) {
    const double r = s->cfg.radius;
    const double W = s->cfg.domain_w;
    const double H = s->cfg.domain_h;
    const double D = s->cfg.domain_d;
    for (int i = 0; i < s->n; i++) {
        if (s->px[i] < r)        { s->px[i] = 2.0 * r       - s->px[i]; s->vx[i] = -s->vx[i]; }
        if (s->px[i] > W - r)    { s->px[i] = 2.0 * (W - r) - s->px[i]; s->vx[i] = -s->vx[i]; }
        if (s->py[i] < r)        { s->py[i] = 2.0 * r       - s->py[i]; s->vy[i] = -s->vy[i]; }
        if (s->py[i] > H - r)    { s->py[i] = 2.0 * (H - r) - s->py[i]; s->vy[i] = -s->vy[i]; }
        if (s->pz[i] < r)        { s->pz[i] = 2.0 * r       - s->pz[i]; s->vz[i] = -s->vz[i]; }
        if (s->pz[i] > D - r)    { s->pz[i] = 2.0 * (D - r) - s->pz[i]; s->vz[i] = -s->vz[i]; }
    }
}

void resolve_pair_collisions(double *vx, double *vy,
                              const double *px, const double *py,
                              const pair_t *pairs, int npairs,
                              double restitution) {
    for (int k = 0; k < npairs; k++) {
        int i = pairs[k].i, j = pairs[k].j;
        double dx = px[j] - px[i];
        double dy = py[j] - py[i];
        double dist2 = dx * dx + dy * dy;
        if (dist2 <= 0.0) continue;                 /* coincident -- normal undefined */
        double inv = 1.0 / sqrt(dist2);
        double nx  = dx * inv;
        double ny  = dy * inv;
        /* v_rel = v_j - v_i; n points i->j; separating iff (v_j - v_i)*n >= 0 */
        double rvx = vx[j] - vx[i];
        double rvy = vy[j] - vy[i];
        double vrel_n = rvx * nx + rvy * ny;
        if (vrel_n >= 0.0) continue;                /* separating; do nothing */
        /* equal masses: J on j is +J*n, on i is -J*n; (1+e)/2 factor */
        double J  = -(1.0 + restitution) * vrel_n * 0.5;
        double Jx = J * nx;
        double Jy = J * ny;
        vx[j] += Jx;  vy[j] += Jy;
        vx[i] -= Jx;  vy[i] -= Jy;
    }
}

void resolve_pair_collisions3d(double *vx, double *vy, double *vz,
                               const double *px, const double *py, const double *pz,
                               const pair_t *pairs, int npairs,
                               double restitution) {
    for (int k = 0; k < npairs; k++) {
        int i = pairs[k].i, j = pairs[k].j;
        double dx = px[j] - px[i];
        double dy = py[j] - py[i];
        double dz = pz[j] - pz[i];
        double dist2 = dx * dx + dy * dy + dz * dz;
        if (dist2 <= 0.0) continue;                 /* coincident -- normal undefined */
        double inv = 1.0 / sqrt(dist2);
        double nx  = dx * inv, ny = dy * inv, nz = dz * inv;
        double rvx = vx[j] - vx[i];
        double rvy = vy[j] - vy[i];
        double rvz = vz[j] - vz[i];
        double vrel_n = rvx * nx + rvy * ny + rvz * nz;
        if (vrel_n >= 0.0) continue;                /* separating; do nothing */
        double J  = -(1.0 + restitution) * vrel_n * 0.5;
        double Jx = J * nx, Jy = J * ny, Jz = J * nz;
        vx[j] += Jx;  vy[j] += Jy;  vz[j] += Jz;
        vx[i] -= Jx;  vy[i] -= Jy;  vz[i] -= Jz;
    }
}

/* Grow s->pairs buffer (used for both candidates and overlapping output) to
 * at least 'needed' entries. Doubles aggressively. */
static void ensure_pair_capacity(sim_t *s, int needed) {
    if (s->pairs_cap >= needed) return;
    while (s->pairs_cap < needed) s->pairs_cap *= 2;
    s->pairs = (pair_t *)realloc(s->pairs, sizeof(pair_t) * (size_t)s->pairs_cap);
    if (!s->pairs) { fprintf(stderr, "sim_step: pair-buffer realloc failed\n"); exit(1); }
}

void sim_step(sim_t *s, sim_metrics_t *out) {
    const double dt  = s->cfg.dt;
    const int    n   = s->n;
    const int    dim = s->cfg.dim;

    struct timespec t0, tg0, tg1, tb1, tn1, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    /* 1a. force pass (semi-implicit Euler).
     *     Dispatches on s->scenario_id; no-op for force-free scenarios. */
    scenario_apply_forces(s);
    /* 1b. position integration with the post-force velocity. */
    if (dim == 3) {
        for (int i = 0; i < n; i++) {
            s->px[i] += s->vx[i] * dt;
            s->py[i] += s->vy[i] * dt;
            s->pz[i] += s->vz[i] * dt;
        }
    } else {
        for (int i = 0; i < n; i++) {
            s->px[i] += s->vx[i] * dt;
            s->py[i] += s->vy[i] * dt;
        }
    }

    /* 2. wall reflection (perfectly elastic) */
    if (dim == 3) wall_reflect3d(s); else wall_reflect(s);

    /* 3. grid build (timed, with internal O(M)/O(N) split for Phase 5) */
    double t_M_in = 0, t_N_in = 0;
    clock_gettime(CLOCK_MONOTONIC, &tg0);
    if (dim == 3)
        grid_build_timed3d(&s->grid, s->px, s->py, s->pz, &t_M_in, &t_N_in);
    else
        grid_build_timed(&s->grid, s->px, s->py, &t_M_in, &t_N_in);
    clock_gettime(CLOCK_MONOTONIC, &tg1);

    /* 4. broad phase: emit candidate pairs (timed) */
    int n_cand;
    while (grid_broad_phase(&s->grid, s->pairs, s->pairs_cap, &n_cand) != 0) {
        ensure_pair_capacity(s, s->pairs_cap * 2);
    }
    clock_gettime(CLOCK_MONOTONIC, &tb1);

    /* Broad-phase outer-iter cost calibration -- separates t_M-prop from
     * t_S-prop in t_broad. Cheap; one walk over cell_start.            */
    long iter_checksum = 0;
    double t_broad_iter = grid_broad_phase_iter_time(&s->grid, &iter_checksum);
    (void)iter_checksum;

    /* 5. narrow phase: circle/sphere overlap filter (in place; timed) */
    int n_overlap;
    if (dim == 3)
        narrow_phase3d(s->pairs, n_cand, s->px, s->py, s->pz, s->cfg.radius,
                       s->pairs, &n_overlap);
    else
        narrow_phase(s->pairs, n_cand, s->px, s->py, s->cfg.radius,
                     s->pairs, &n_overlap);
    clock_gettime(CLOCK_MONOTONIC, &tn1);

    /* 6. resolve overlapping-pair collisions */
    if (dim == 3)
        resolve_pair_collisions3d(s->vx, s->vy, s->vz, s->px, s->py, s->pz,
                                  s->pairs, n_overlap, s->cfg.restitution);
    else
        resolve_pair_collisions(s->vx, s->vy, s->px, s->py,
                                s->pairs, n_overlap, s->cfg.restitution);

    clock_gettime(CLOCK_MONOTONIC, &t1);

    s->last_cand_count = n_cand;
    s->last_pair_count = n_overlap;
    s->frame_index++;

    if (out) {
        out->cell_size        = s->cfg.cell_size;
        out->N                = n;
        out->num_cells        = s->grid.n_cells;
        out->S                = grid_S(&s->grid);
        out->candidate_pairs  = n_cand;
        out->collisions       = n_overlap;
        out->t_grid           = ts_diff_s(&tg0, &tg1);
        out->t_broad          = ts_diff_s(&tg1, &tb1);
        out->t_narrow         = ts_diff_s(&tb1, &tn1);
        out->t_total          = ts_diff_s(&t0,  &t1);
        out->t_M              = t_M_in;
        out->t_N              = t_N_in;
        out->t_broad_iter     = t_broad_iter;
        out->kinetic_energy   = sim_total_kinetic_energy(s);
        if (dim == 3) {
            double mx = 0, my = 0, mz = 0;
            for (int i = 0; i < n; i++) { mx += s->vx[i]; my += s->vy[i]; mz += s->vz[i]; }
            out->momentum_magnitude = sqrt(mx * mx + my * my + mz * mz);
        } else {
            double mpx, mpy;
            sim_total_momentum(s, &mpx, &mpy);
            out->momentum_magnitude = sqrt(mpx * mpx + mpy * mpy);
        }
    }
}

double sim_total_kinetic_energy(const sim_t *s) {
    /* unit masses */
    double ke = 0.0;
    if (s->cfg.dim == 3) {
        for (int i = 0; i < s->n; i++) {
            ke += 0.5 * (s->vx[i] * s->vx[i] + s->vy[i] * s->vy[i] + s->vz[i] * s->vz[i]);
        }
    } else {
        for (int i = 0; i < s->n; i++) {
            ke += 0.5 * (s->vx[i] * s->vx[i] + s->vy[i] * s->vy[i]);
        }
    }
    return ke;
}

void sim_total_momentum(const sim_t *s, double *px_total, double *py_total) {
    double pxt = 0.0, pyt = 0.0;
    for (int i = 0; i < s->n; i++) { pxt += s->vx[i]; pyt += s->vy[i]; }
    *px_total = pxt;
    *py_total = pyt;
}
