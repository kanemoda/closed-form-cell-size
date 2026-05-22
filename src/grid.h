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
    int    gw, gh;       /* grid width / height in cells */
    int    n_cells;      /* = gw * gh                     */
    int    n;            /* particle count                */
    int   *cell_count;   /* size n_cells                  */
    int   *cell_start;   /* size n_cells + 1              */
    int   *sorted_idx;   /* size n: particle ids by cell  */
    int   *cell_idx;     /* size n: cell id per particle  */
} grid_t;

void grid_init (grid_t *g, int n, double domain_w, double domain_h, double cell_size);
void grid_free (grid_t *g);
void grid_build(grid_t *g, const double *px, const double *py);

/* M(l) = number of cells, S(l) = sum of n_i^2 over cells. */
int  grid_M(const grid_t *g);
long grid_S(const grid_t *g);

/*
 * Broad phase: emit ALL candidate pairs from intra-cell + half-stencil,
 * with pair indices normalised so i < j. No distance check.
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

#endif
