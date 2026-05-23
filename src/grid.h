#ifndef GRID_H
#define GRID_H

#include <stddef.h>

typedef struct { int i, j; } pair_t;   /* always normalised so i < j */

/*
 * Uniform-grid broad phase. Counting-sort binning:
 *   1. compute cell index per particle
 *   2. per-cell counts
 *   3. exclusive prefix sum -> cell_start
 *   4. scatter into sorted_idx
 *
 * Query (broad phase): per cell, emit intra-cell pairs + pairs with
 * particles in the four half-neighbor cells (offsets (+1,0), (-1,+1),
 * (0,+1), (+1,+1) in 2D). No distance test -- that is the narrow phase.
 *
 * Correctness: for the 1-ring half-stencil to be complete (no missed
 * overlapping pairs), cell_size must be >= 2*radius.
 */
typedef struct {
    double cell_size;
    int    dim;          /* 2 or 3 (set by grid_init / grid_init3d)  */
    int    gw, gh, gd;   /* grid extent in cells (gd = 1 when dim==2)*/
    int    n_cells;      /* = gw * gh * gd                           */
    int    n;            /* particle count                */
    int   *cell_count;   /* size n_cells                  */
    int   *cell_start;   /* size n_cells + 1              */
    int   *sorted_idx;   /* size n: particle ids by cell  */
    int   *cell_idx;     /* size n: cell id per particle  */
    long   S_cached;     /* Sum n_i^2 from the last build, computed
                            incrementally during the count pass     */
    double domain_w;     /* kept for grid_resize / grid_count_S_at  */
    double domain_h;
    double domain_d;     /* z extent (3D only; 0 when dim==2)        */
} grid_t;

void grid_init (grid_t *g, int n, double domain_w, double domain_h, double cell_size);

/* 3D counterpart: cubic-ish cells of edge cell_size over a domain_w x
 * domain_h x domain_d box. n_cells = gw*gh*gd. Sets dim=3. */
void grid_init3d(grid_t *g, int n, double domain_w, double domain_h,
                 double domain_d, double cell_size);
void grid_free (grid_t *g);
void grid_build(grid_t *g, const double *px, const double *py);

/* Like grid_build, but fills *t_M (the strictly O(M) sub-pass: cell-count
 * memset + prefix sum + offset memset) and *t_N (the strictly O(N)
 * sub-passes: cell_idx compute + count+S accumulation + scatter). Either
 * pointer may be NULL. Used by the adaptive controller, which wants
 * t_M / M and t_S / S as the per-cell and per-pair coefficient inputs to
 * its EMA. */
void grid_build_timed(grid_t *g, const double *px, const double *py,
                      double *t_M, double *t_N);

/* 3D build (needs the z coordinate). Same counting-sort pipeline; only the
 * per-particle cell index differs (cz*gh*gw + cy*gw + cx). */
void grid_build3d(grid_t *g, const double *px, const double *py, const double *pz);
void grid_build_timed3d(grid_t *g, const double *px, const double *py,
                        const double *pz, double *t_M, double *t_N);

/* Resize the grid in place to a new cell size (reallocates cell arrays as
 * needed). After grid_resize you must call grid_build again. Used by the
 * adaptive controller; cheap when n_cells does not change. */
void grid_resize(grid_t *g, double new_cell_size);

/* M(l) = number of cells, S(l) = sum of n_i^2 over cells.
 * grid_S returns the value cached during the last grid_build (O(1));
 * no extra pass over cells. */
int  grid_M(const grid_t *g);
long grid_S(const grid_t *g);

/* Count-only occupancy probe at an arbitrary cell size: runs only the
 * count pass for a fresh grid at 'ell' over the same particles and
 * returns S = sum n_i^2. No prefix sum, no scatter, no sorted_idx;
 * O(N + M(ell)) but the M term is just one memset over a scratch
 * counts array. Used by the Mode B D2 probe at ell/gamma and ell*gamma.
 *
 * Allocates internally; safe with concurrent grid state (does not touch g). */
long grid_count_S_at(double domain_w, double domain_h, int n,
                     const double *px, const double *py, double ell);

/* 3D count-only S probe (Mode B D2 estimator in 3D). */
long grid_count_S_at3d(double domain_w, double domain_h, double domain_d, int n,
                       const double *px, const double *py, const double *pz,
                       double ell);

/*
 * Broad phase: emit ALL candidate pairs from intra-cell + half-stencil,
 * with pair indices normalised so i < j. No distance check. Dispatches on
 * g->dim: 4-neighbor half-stencil in 2D, 13-neighbor half-stencil in 3D.
 *
 * Returns 0 on success, -1 if pairs_out has insufficient capacity (in
 * which case *count_out is the number of pairs emitted before overflow).
 */
int grid_broad_phase(const grid_t *g,
                      pair_t *pairs_out, int max_pairs, int *count_out);

/*
 * Narrow phase: from a buffer of candidate pairs, keep those whose two
 * particles actually overlap (distance < 2r). Safe to use with out==in
 * (compacts in place).
 */
void narrow_phase(const pair_t *cands, int n_cands,
                   const double *px, const double *py, double r,
                   pair_t *out, int *n_overlap);

/* 3D sphere-sphere overlap filter (distance < 2r in 3D). */
void narrow_phase3d(const pair_t *cands, int n_cands,
                    const double *px, const double *py, const double *pz, double r,
                    pair_t *out, int *n_overlap);

/*
 * Convenience: broad + narrow in one call. Same buffer used as candidate
 * scratch and final output. Returns 0 on success, -1 on capacity overflow.
 */
int grid_find_overlapping_pairs(const grid_t *g,
                                 const double *px, const double *py,
                                 double r,
                                 pair_t *pairs_out, int max_pairs, int *count_out);

/* Same overlap set, but via O(N^2) brute force. Reference oracle. */
int bruteforce_find_overlapping_pairs(int n,
                                       const double *px, const double *py,
                                       double r,
                                       pair_t *pairs_out, int max_pairs, int *count_out);

/* 3D convenience + brute-force reference (used by the Phase 6 correctness gate). */
int grid_find_overlapping_pairs3d(const grid_t *g,
                                   const double *px, const double *py, const double *pz,
                                   double r,
                                   pair_t *pairs_out, int max_pairs, int *count_out);
int bruteforce_find_overlapping_pairs3d(int n,
                                        const double *px, const double *py, const double *pz,
                                        double r,
                                        pair_t *pairs_out, int max_pairs, int *count_out);

/* Phase 5: time the broad-phase OUTER CELL ITERATION ONLY, with the same
 * memory-access pattern as grid_broad_phase but no pair emission. Used to
 * split t_broad into its O(M) component (which belongs in the controller's
 * t_M) and its O(S) component (the inner pair-emit, which belongs in t_S).
 * Returns the elapsed wall-clock (s) of the dryrun. The integer return
 * value (a checksum) is opaque -- only there to defeat dead-code elimination. */
double grid_broad_phase_iter_time(const grid_t *g, long *checksum_out);

#endif
