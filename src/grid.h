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
 * Query: per cell, intra-cell pairs + half-neighbor cells
 * (offsets (+1,0), (-1,+1), (0,+1), (+1,+1) in 2D).
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
 * Fill pairs_out with all *truly overlapping* pairs (broad + narrow phase).
 * Returns 0 on success, -1 if pairs_out has insufficient capacity
 * (in which case *count_out is the number of pairs found before overflow).
 */
int grid_find_overlapping_pairs(const grid_t *g,
                                 const double *px, const double *py,
                                 double r,
                                 pair_t *pairs_out, int max_pairs, int *count_out);

/* Same, but via O(N^2) brute force. Reference oracle. */
int bruteforce_find_overlapping_pairs(int n,
                                       const double *px, const double *py,
                                       double r,
                                       pair_t *pairs_out, int max_pairs, int *count_out);

#endif
