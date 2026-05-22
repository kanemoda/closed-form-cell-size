#ifndef SIM_H
#define SIM_H

#include "rng.h"
#include "config.h"
#include "grid.h"

/*
 * Single-threaded 2D particle simulator (Phase 1).
 *
 * Per step:
 *   1. positions <- positions + velocities * dt   (explicit Euler)
 *   2. wall reflection (perfectly elastic)
 *   3. grid build  (counting sort)
 *   4. find overlapping pairs (broad + narrow phase)
 *   5. resolve pair collisions (velocity-only, elastic with restitution e)
 *
 * Determinism: with a fixed config + seed, every internal step is a
 * deterministic function of inputs -- iteration orders are fixed, no
 * floating-point ops depend on memory layout. Same seed -> same bits.
 */
typedef struct {
    config_t     cfg;
    int          n;
    double      *px, *py;
    double      *vx, *vy;
    rng_state_t  rng;
    grid_t       grid;
    pair_t      *pairs;
    int          pairs_cap;
    int          last_pair_count;
} sim_t;

void sim_init       (sim_t *s, const config_t *cfg);
void sim_init_random(sim_t *s);
void sim_free       (sim_t *s);
void sim_step       (sim_t *s);

double sim_total_kinetic_energy(const sim_t *s);
void   sim_total_momentum      (const sim_t *s, double *px_total, double *py_total);

/* Apply elastic-with-restitution velocity impulses to each overlapping pair.
 * Equal masses; n points from i to j; impulse skipped if pair is separating. */
void resolve_pair_collisions(double *vx, double *vy,
                              const double *px, const double *py,
                              const pair_t *pairs, int npairs,
                              double restitution);

#endif
