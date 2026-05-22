#include "grid.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

void grid_init(grid_t *g, int n, double domain_w, double domain_h, double cell_size) {
    g->cell_size  = cell_size;
    g->gw         = (int)ceil(domain_w / cell_size);
    g->gh         = (int)ceil(domain_h / cell_size);
    if (g->gw < 1) g->gw = 1;
    if (g->gh < 1) g->gh = 1;
    g->n_cells    = g->gw * g->gh;
    g->n          = n;
    g->cell_count = (int *)calloc((size_t)g->n_cells,     sizeof(int));
    g->cell_start = (int *)calloc((size_t)g->n_cells + 1, sizeof(int));
    g->sorted_idx = (int *)calloc((size_t)n,              sizeof(int));
    g->cell_idx   = (int *)calloc((size_t)n,              sizeof(int));
}

void grid_free(grid_t *g) {
    free(g->cell_count);
    free(g->cell_start);
    free(g->sorted_idx);
    free(g->cell_idx);
    g->cell_count = g->cell_start = g->sorted_idx = g->cell_idx = NULL;
}

void grid_build(grid_t *g, const double *px, const double *py) {
    const int    gw  = g->gw, gh = g->gh;
    const double inv = 1.0 / g->cell_size;
    const int    n   = g->n;

    /* 1. cell index per particle (clamped) */
    for (int i = 0; i < n; i++) {
        int cx = (int)floor(px[i] * inv);
        int cy = (int)floor(py[i] * inv);
        if (cx < 0)   cx = 0;
        if (cx >= gw) cx = gw - 1;
        if (cy < 0)   cy = 0;
        if (cy >= gh) cy = gh - 1;
        g->cell_idx[i] = cy * gw + cx;
    }

    /* 2. counts */
    memset(g->cell_count, 0, sizeof(int) * (size_t)g->n_cells);
    for (int i = 0; i < n; i++) g->cell_count[g->cell_idx[i]]++;

    /* 3. exclusive prefix sum, with sentinel cell_start[n_cells] = n */
    int run = 0;
    for (int c = 0; c < g->n_cells; c++) {
        g->cell_start[c] = run;
        run += g->cell_count[c];
    }
    g->cell_start[g->n_cells] = run;

    /* 4. scatter -- reuse cell_count as a per-cell insertion offset.
     *    After the loop, offset[c] equals the original cell_count[c],
     *    so cell_count is naturally restored.                       */
    int *offset = g->cell_count;
    memset(offset, 0, sizeof(int) * (size_t)g->n_cells);
    for (int i = 0; i < n; i++) {
        int c = g->cell_idx[i];
        g->sorted_idx[g->cell_start[c] + offset[c]] = i;
        offset[c]++;
    }
}

int  grid_M(const grid_t *g) { return g->n_cells; }

long grid_S(const grid_t *g) {
    long s = 0;
    for (int c = 0; c < g->n_cells; c++) {
        long ni = g->cell_count[c];
        s += ni * ni;
    }
    return s;
}

/* test pair (ia, ib): if overlapping, add to out in normalised (i<j) order. */
static inline int try_add_pair(const double *px, const double *py, double rr,
                                int ia, int ib,
                                pair_t *out, int max_pairs, int *count) {
    double dx = px[ia] - px[ib];
    double dy = py[ia] - py[ib];
    if (dx * dx + dy * dy < rr) {
        if (*count >= max_pairs) return -1;
        int i = ia < ib ? ia : ib;
        int j = ia < ib ? ib : ia;
        out[*count].i = i;
        out[*count].j = j;
        (*count)++;
    }
    return 0;
}

int grid_find_overlapping_pairs(const grid_t *g,
                                 const double *px, const double *py,
                                 double r,
                                 pair_t *pairs_out, int max_pairs, int *count_out) {
    *count_out = 0;
    const double rr = (2.0 * r) * (2.0 * r);
    const int    gw = g->gw, gh = g->gh;

    /* Half-stencil: (+1,0), (-1,+1), (0,+1), (+1,+1).
     * Each unordered pair of adjacent cells is visited from exactly one side. */
    static const int dxs[4] = {  1, -1,  0,  1 };
    static const int dys[4] = {  0,  1,  1,  1 };

    for (int cy = 0; cy < gh; cy++) {
        for (int cx = 0; cx < gw; cx++) {
            int c  = cy * gw + cx;
            int a0 = g->cell_start[c];
            int a1 = g->cell_start[c + 1];

            /* intra-cell */
            for (int a = a0; a < a1; a++) {
                int ia = g->sorted_idx[a];
                for (int b = a + 1; b < a1; b++) {
                    int ib = g->sorted_idx[b];
                    if (try_add_pair(px, py, rr, ia, ib, pairs_out, max_pairs, count_out) != 0)
                        return -1;
                }
            }

            /* half-neighbor cells */
            for (int k = 0; k < 4; k++) {
                int nx = cx + dxs[k];
                int ny = cy + dys[k];
                if (nx < 0 || nx >= gw || ny < 0 || ny >= gh) continue;
                int nc = ny * gw + nx;
                int b0 = g->cell_start[nc];
                int b1 = g->cell_start[nc + 1];
                for (int a = a0; a < a1; a++) {
                    int ia = g->sorted_idx[a];
                    for (int b = b0; b < b1; b++) {
                        int ib = g->sorted_idx[b];
                        if (try_add_pair(px, py, rr, ia, ib, pairs_out, max_pairs, count_out) != 0)
                            return -1;
                    }
                }
            }
        }
    }
    return 0;
}

int bruteforce_find_overlapping_pairs(int n,
                                       const double *px, const double *py,
                                       double r,
                                       pair_t *pairs_out, int max_pairs, int *count_out) {
    *count_out = 0;
    const double rr = (2.0 * r) * (2.0 * r);
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double dx = px[j] - px[i];
            double dy = py[j] - py[i];
            if (dx * dx + dy * dy < rr) {
                if (*count_out >= max_pairs) return -1;
                pairs_out[*count_out].i = i;
                pairs_out[*count_out].j = j;
                (*count_out)++;
            }
        }
    }
    return 0;
}
