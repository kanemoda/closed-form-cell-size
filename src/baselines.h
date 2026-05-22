#ifndef BASELINES_H
#define BASELINES_H

#include "config.h"

/*
 * Fixed-ℓ baselines and the static-oracle driver (spec §9.2).
 *
 * Two fixed-ℓ baselines are provided:
 *   - Ericson:    ℓ = 2r  (the half-stencil correctness floor)
 *   - Analytical: ℓ = (V/N)^(1/d), clamped to ≥ 2r
 *     Source: Frenkel & Smit, *Understanding Molecular Simulation* (2002),
 *     §5.2.2 -- the standard cell-list balance argument under uniform
 *     density with equal cell- and pair-cost coefficients. The spec asks
 *     for "Sigurgeirsson analytical formula" but leaves the exact formula
 *     unspecified (see TODO.md item 2); we use the F&S formula as the
 *     standard published heuristic until the Sigurgeirsson reference is
 *     resolved.
 *
 * The static oracle runs the full simulation once per candidate ℓ
 * (deterministic given seed+config) and reports the ℓ that minimises
 * mean detection cost (mean of t_grid + t_broad + t_narrow over frames).
 */

static inline double baseline_ericson(double radius) { return 2.0 * radius; }

double baseline_uniform_analytical(int N, double domain_w, double domain_h,
                                   double radius);

/* Aggregated stats from one full scenario run. */
typedef struct {
    double mean_t_grid;
    double mean_t_broad;
    double mean_t_narrow;
    double mean_t_detect;
    long   total_candidates;
    long   total_collisions;
    int    frames;
} sim_run_summary_t;

/* Run the configured scenario with cell_size forced to the given value.
 * Disables CSV logging during the run. */
void run_sim_with_cell_size(const config_t *cfg, double cell_size,
                            sim_run_summary_t *out);

/* Static oracle: try each ell in 'ells' (calling run_sim_with_cell_size
 * for each), keep the one with minimum mean_t_detect. */
typedef struct {
    int     n_candidates;
    double *ells;          /* copied from input                       */
    double *mean_detects;  /* mean_t_detect per ell                   */
    int     best_index;
    double  best_ell;
    double  best_mean_detect;
} static_oracle_result_t;

void static_oracle_run (const config_t *cfg, const double *ells, int n_ells,
                        static_oracle_result_t *out);
void static_oracle_free(static_oracle_result_t *r);

#endif
