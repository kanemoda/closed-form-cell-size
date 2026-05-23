#ifndef ORACLE_H
#define ORACLE_H

/*
 * Phase 4 oracle harness.
 *
 * Given a *frozen* single-frame particle snapshot (px, py), evaluate the
 * true broad-phase detection cost at a configurable geometric sweep of
 * candidate cell sizes ℓ. For each ℓ, build a fresh grid + broad + narrow
 * (matching sim_step's pipeline) and measure t_grid, t_broad, t_narrow via
 * CLOCK_MONOTONIC; takes the median over K timing repeats to suppress
 * wall-clock noise. Also records S, num_cells, candidate_pairs at that ℓ.
 *
 * The sweep is non-perturbing: it never touches the caller's sim state.
 *
 * Used by:
 *   - per-frame oracle (best ℓ each frame -- the adaptive upper bound),
 *     invoked from main.c after each sim_step;
 *   - static oracle (best single fixed ℓ over a run), in baselines.c.
 */

typedef struct {
    double ell;
    int    num_cells;        /* M(ℓ)                     */
    long   S;                /* Σ n_i²                   */
    int    candidate_pairs;  /* broad-phase output count */
    double t_grid;           /* median over K repeats    */
    double t_broad;          /* median over K repeats    */
    double t_narrow;         /* median over K repeats    */
    double t_detect;         /* t_grid + t_broad + t_narrow */
} oracle_point_t;

/* Generate ells[k] = ell_min * gamma^k for k in [0, n), stopping at
 * ell > ell_max OR n == max_pts. Returns n. */
int  oracle_make_geo_sweep(double ell_min, double ell_max, double gamma,
                           double *ells, int max_pts);

/* Evaluate one (ell, snapshot). Allocates internal scratch.
 * K = number of timing repeats; the median is taken per phase. */
void oracle_eval(int n, const double *px, const double *py, double radius,
                 double domain_w, double domain_h,
                 double ell, int K,
                 oracle_point_t *out);

/* Evaluate the full sweep on a snapshot. Caller-allocated 'out' of size n_ells. */
void oracle_eval_sweep(int n, const double *px, const double *py, double radius,
                       double domain_w, double domain_h,
                       const double *ells, int n_ells, int K,
                       oracle_point_t *out);

/* Split-sample evaluation (de-biasing). Times the full detect pipeline
 * (build+broad+narrow) 2K times on the SAME (snapshot, ell), interleaved into
 * two independent groups, and returns each group's median as two independent
 * t_detect estimates. The per-frame oracle SELECTS argmin over 'out_sel' and
 * REPORTS the corresponding 'out_eval', so the reported cost is not the
 * optimistically-biased minimum of the same noisy sample it was chosen on.
 * out_sel / out_eval share ell, S, num_cells; only the timings differ. */
void oracle_eval_split(int n, const double *px, const double *py, double radius,
                       double domain_w, double domain_h,
                       double ell, int K,
                       oracle_point_t *out_sel, oracle_point_t *out_eval);

/* 3D counterparts (3D grid, sphere narrow phase). Same semantics as the 2D
 * versions; the split-sample variant feeds the 3D comparison's oracles. */
void oracle_eval3d(int n, const double *px, const double *py, const double *pz,
                   double radius, double domain_w, double domain_h, double domain_d,
                   double ell, int K, oracle_point_t *out);
void oracle_eval_sweep3d(int n, const double *px, const double *py, const double *pz,
                         double radius, double domain_w, double domain_h, double domain_d,
                         const double *ells, int n_ells, int K, oracle_point_t *out);
void oracle_eval_split3d(int n, const double *px, const double *py, const double *pz,
                         double radius, double domain_w, double domain_h, double domain_d,
                         double ell, int K,
                         oracle_point_t *out_sel, oracle_point_t *out_eval);

#endif
