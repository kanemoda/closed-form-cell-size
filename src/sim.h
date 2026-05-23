#ifndef SIM_H
#define SIM_H

#include "rng.h"
#include "config.h"
#include "grid.h"

/*
 * Single-threaded 2D particle simulator.
 *
 * Integrator: semi-implicit (symplectic) Euler:
 *   1. v <- v + (F/m) * dt          [no forces in Phase 1; no-op]
 *   2. x <- x + v_new * dt
 * With F=0 this is bit-identical to explicit Euler (Phase 1), but the
 * structure is in place so Phase 3 can add gravity/springs and keep
 * total energy bounded.
 *
 * Per step (timed regions marked *):
 *   1. integrator (force-pass then position update)
 *   2. wall reflection (perfectly elastic)
 * * 3. grid build   (counting sort)              -> t_grid
 * * 4. broad phase  (emit candidate pairs)        -> t_broad
 * * 5. narrow phase (circle-circle overlap test)  -> t_narrow
 *   6. resolve pair collisions (elastic, restitution e)
 *   Whole frame elapsed -> t_total.
 */

/* Per-frame metrics. Populated by sim_step when 'out' is non-NULL.
 * 'frame' is *not* touched by sim_step -- the caller sets it. */
typedef struct {
    int    frame;
    double cell_size;
    int    N;
    int    num_cells;             /* grid.n_cells (model: M = V/ell^d) */
    long   S;                     /* sum of n_i^2 over cells (Renyi-2 building block) */
    int    candidate_pairs;       /* broad-phase output count */
    int    collisions;            /* narrow-phase survivors / impulse count */
    double t_grid;                /* seconds; total grid_build time = t_M + t_N */
    double t_broad;
    double t_narrow;
    double t_total;
    /* Phase 5: spec-mandated split of grid build into O(M) and O(N) work.
     *   t_M = grid_build's strictly-O(M) sub-passes (cell-count clear,
     *         prefix sum, offset clear). Per-cell coefficient signal.
     *   t_N = grid_build's strictly-O(N) sub-passes (cell-index compute,
     *         count + S accumulation, scatter). Insertion cost; logged
     *         but not used by the closed-form formula (which only cares
     *         about ell-dependent terms).
     *   The pair-test coefficient signal for the controller is
     *   t_broad + t_narrow (spec calls this t_S). */
    double t_M;
    double t_N;
    /* Broad-phase outer-loop overhead (O(M)), measured by a dryrun pass
     * over cell_start without pair emission. The controller adds this to
     * t_M (because it scales with M, not S) and subtracts it from t_broad
     * to recover the actual pair-emit cost. */
    double t_broad_iter;
    double kinetic_energy;
    double momentum_magnitude;
} sim_metrics_t;

typedef struct {
    config_t     cfg;
    int          n;
    double      *px, *py;
    double      *vx, *vy;
    double      *pz, *vz;          /* 3D only (NULL when cfg.dim == 2)      */
    rng_state_t  rng;
    grid_t       grid;
    pair_t      *pairs;            /* doubles as broad candidate buffer + narrow output */
    int          pairs_cap;
    int          last_cand_count;
    int          last_pair_count;

    /* Phase 3: scenario state. Set by scenario_init(). The force-field
     * dispatcher (scenario_apply_forces) reads scenario_id; sim_step()
     * increments frame_index at the end of each step so time-dependent
     * forces can use t = frame_index * dt. */
    int          scenario_id;      /* index into the scenario vtable; default 0 */
    int          frame_index;      /* steps elapsed since scenario_init        */
} sim_t;

void sim_init       (sim_t *s, const config_t *cfg);
void sim_init_random(sim_t *s);
void sim_free       (sim_t *s);

/* One integrator step. If 'out' is non-NULL, per-frame metrics are filled.
 * The caller is responsible for setting out->frame; sim_step touches all
 * other fields. */
void sim_step       (sim_t *s, sim_metrics_t *out);

double sim_total_kinetic_energy(const sim_t *s);
void   sim_total_momentum      (const sim_t *s, double *px_total, double *py_total);

/* Apply elastic-with-restitution velocity impulses to each overlapping pair.
 * Equal masses; n points from i to j; impulse skipped if pair is separating. */
void resolve_pair_collisions(double *vx, double *vy,
                              const double *px, const double *py,
                              const pair_t *pairs, int npairs,
                              double restitution);

/* 3D elastic-with-restitution impulse resolution (equal masses). */
void resolve_pair_collisions3d(double *vx, double *vy, double *vz,
                               const double *px, const double *py, const double *pz,
                               const pair_t *pairs, int npairs,
                               double restitution);

#endif
