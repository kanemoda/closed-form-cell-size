/*
 * Verification gate 1 (spec section 10):
 *   "The set of overlapping pairs detected through the grid pipeline
 *    exactly equals the set from an O(N^2) brute-force check, on small N
 *    -- no missed collisions, no spurious ones."
 *
 * We test several N, several cell sizes (always >= 2r), non-overlapping
 * AND overlapping placements, and the boundary case cell_size == 2r exactly.
 */

#include "rng.h"
#include "grid.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int cmp_pair(const void *a, const void *b) {
    const pair_t *pa = (const pair_t *)a;
    const pair_t *pb = (const pair_t *)b;
    if (pa->i != pb->i) return pa->i - pb->i;
    return pa->j - pb->j;
}

static void place_random_no_overlap(rng_state_t *rng,
                                     int n, double W, double H, double r,
                                     double *px, double *py) {
    double rr_min = (2 * r) * (2 * r);
    for (int i = 0; i < n; i++) {
        int placed = 0;
        for (int t = 0; t < 20000; t++) {
            double x = rng_uniform(rng, r, W - r);
            double y = rng_uniform(rng, r, H - r);
            int ok = 1;
            for (int j = 0; j < i; j++) {
                double dx = x - px[j], dy = y - py[j];
                if (dx * dx + dy * dy < rr_min) { ok = 0; break; }
            }
            if (ok) { px[i] = x; py[i] = y; placed = 1; break; }
        }
        if (!placed) {
            fprintf(stderr, "  place_random_no_overlap: failed at i=%d (density too high)\n", i);
            exit(1);
        }
    }
}

static void place_random_allow_overlap(rng_state_t *rng,
                                        int n, double W, double H, double r,
                                        double *px, double *py) {
    for (int i = 0; i < n; i++) {
        px[i] = rng_uniform(rng, r, W - r);
        py[i] = rng_uniform(rng, r, H - r);
    }
}

static int run_case(int N, double W, double H, double r,
                    double cell_size, uint64_t seed, int allow_overlap) {
    rng_state_t rng;
    rng_seed(&rng, seed);

    double *px = (double *)malloc(sizeof(double) * (size_t)N);
    double *py = (double *)malloc(sizeof(double) * (size_t)N);

    if (allow_overlap) place_random_allow_overlap(&rng, N, W, H, r, px, py);
    else                place_random_no_overlap   (&rng, N, W, H, r, px, py);

    grid_t g;
    grid_init(&g, N, W, H, cell_size);
    grid_build(&g, px, py);

    int max_pairs = N * (N - 1) / 2 + 1;
    pair_t *p_grid  = (pair_t *)malloc(sizeof(pair_t) * (size_t)max_pairs);
    pair_t *p_brute = (pair_t *)malloc(sizeof(pair_t) * (size_t)max_pairs);
    int n_grid = 0, n_brute = 0;

    int rg = grid_find_overlapping_pairs   (&g, px, py, r, p_grid,  max_pairs, &n_grid);
    int rb = bruteforce_find_overlapping_pairs(N, px, py, r, p_brute, max_pairs, &n_brute);
    if (rg != 0 || rb != 0) {
        fprintf(stderr, "  case N=%d cell=%g: BUFFER OVERFLOW (shouldn't happen)\n", N, cell_size);
        free(p_grid); free(p_brute); grid_free(&g); free(px); free(py);
        return 0;
    }

    qsort(p_grid,  (size_t)n_grid,  sizeof(pair_t), cmp_pair);
    qsort(p_brute, (size_t)n_brute, sizeof(pair_t), cmp_pair);

    int ok = 1;
    if (n_grid != n_brute) {
        fprintf(stderr, "  MISMATCH counts: grid=%d brute=%d\n", n_grid, n_brute);
        ok = 0;
    } else {
        for (int k = 0; k < n_grid; k++) {
            if (p_grid[k].i != p_brute[k].i || p_grid[k].j != p_brute[k].j) {
                fprintf(stderr, "  MISMATCH @ %d: grid=(%d,%d) brute=(%d,%d)\n",
                        k, p_grid[k].i, p_grid[k].j, p_brute[k].i, p_brute[k].j);
                ok = 0;
                break;
            }
        }
    }

    fprintf(stderr,
        "  N=%4d  W=%5g H=%5g r=%g cell=%4.2f seed=%2llu overlap=%d pairs=%4d  %s\n",
        N, W, H, r, cell_size,
        (unsigned long long)seed, allow_overlap, n_brute,
        ok ? "OK" : "FAIL");

    free(p_grid); free(p_brute);
    grid_free(&g);
    free(px); free(py);
    return ok;
}

int main(void) {
    fprintf(stderr, "=== Gate 1: grid overlapping-pair set == brute-force set ===\n");
    int all_ok = 1;

    /* Non-overlapping placements: pair count is 0; checks that
     * the grid pipeline doesn't manufacture spurious pairs. */
    all_ok &= run_case(  50, 30,  30, 0.5, 1.0,  1, 0);
    all_ok &= run_case( 100, 50,  50, 0.5, 1.5,  2, 0);
    all_ok &= run_case( 200, 80,  80, 0.5, 2.0,  3, 0);
    all_ok &= run_case( 200, 80,  80, 0.5, 4.0,  4, 0);
    all_ok &= run_case( 500,100, 100, 0.5, 1.0,  5, 0);
    all_ok &= run_case( 500,100, 100, 0.5, 8.0,  6, 0);

    /* Dense placements (allow overlap): pair count is non-zero;
     * checks that the grid pipeline doesn't miss any. */
    all_ok &= run_case( 100,  10,  10, 0.5, 1.0, 10, 1);   /* cell == 2r boundary */
    all_ok &= run_case( 100,  10,  10, 0.5, 1.5, 11, 1);
    all_ok &= run_case( 100,  10,  10, 0.5, 2.0, 12, 1);
    all_ok &= run_case( 200,  15,  15, 0.5, 1.0, 13, 1);
    all_ok &= run_case( 200,  15,  15, 0.5, 2.5, 14, 1);
    all_ok &= run_case( 500,  20,  20, 0.5, 1.0, 15, 1);
    all_ok &= run_case( 500,  20,  20, 0.5, 4.0, 16, 1);
    all_ok &= run_case(1000,  30,  30, 0.5, 1.0, 17, 1);
    all_ok &= run_case(1000,  30,  30, 0.5, 6.0, 18, 1);

    /* Non-square / non-uniform domain shapes. */
    all_ok &= run_case( 300, 12,  60, 0.5, 1.0, 20, 1);
    all_ok &= run_case( 300, 60,  12, 0.5, 1.0, 21, 1);
    all_ok &= run_case( 300, 12,  60, 0.5, 3.0, 22, 1);

    /* Vary radius. */
    all_ok &= run_case( 300,  20,  20, 0.3, 0.6, 30, 1);
    all_ok &= run_case( 300,  20,  20, 0.7, 1.4, 31, 1);
    all_ok &= run_case( 300,  20,  20, 0.7, 5.0, 32, 1);

    fprintf(stderr, "\n%s\n", all_ok ? "Gate 1: ALL OK" : "Gate 1: FAILED");
    return all_ok ? 0 : 1;
}
