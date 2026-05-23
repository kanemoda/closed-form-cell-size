#define _POSIX_C_SOURCE 199309L
#include "grid.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

static inline double ts_diff_s(const struct timespec *a, const struct timespec *b) {
    return (double)(b->tv_sec - a->tv_sec) + (double)(b->tv_nsec - a->tv_nsec) * 1e-9;
}

void grid_init(grid_t *g, int n, double domain_w, double domain_h, double cell_size) {
    g->cell_size  = cell_size;
    g->domain_w   = domain_w;
    g->domain_h   = domain_h;
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
    g->S_cached   = 0;
}

void grid_resize(grid_t *g, double new_cell_size) {
    if (new_cell_size == g->cell_size) return;
    g->cell_size = new_cell_size;
    int new_gw = (int)ceil(g->domain_w / new_cell_size);
    int new_gh = (int)ceil(g->domain_h / new_cell_size);
    if (new_gw < 1) new_gw = 1;
    if (new_gh < 1) new_gh = 1;
    int new_cells = new_gw * new_gh;
    if (new_cells != g->n_cells) {
        free(g->cell_count); free(g->cell_start);
        g->cell_count = (int *)calloc((size_t)new_cells,     sizeof(int));
        g->cell_start = (int *)calloc((size_t)new_cells + 1, sizeof(int));
    }
    g->gw = new_gw; g->gh = new_gh; g->n_cells = new_cells;
    /* sorted_idx and cell_idx are sized by n, no need to realloc. */
    /* S_cached is stale until the next grid_build. */
    g->S_cached = 0;
}

void grid_free(grid_t *g) {
    free(g->cell_count);
    free(g->cell_start);
    free(g->sorted_idx);
    free(g->cell_idx);
    g->cell_count = g->cell_start = g->sorted_idx = g->cell_idx = NULL;
}

void grid_build(grid_t *g, const double *px, const double *py) {
    grid_build_timed(g, px, py, NULL, NULL);
}

void grid_build_timed(grid_t *g, const double *px, const double *py,
                      double *t_M_out, double *t_N_out) {
    const int    gw  = g->gw, gh = g->gh;
    const double inv = 1.0 / g->cell_size;
    const int    n   = g->n;

    /* The build interleaves strictly-O(N) and strictly-O(M) sub-passes;
     * we accumulate each into its own bucket so the adapter sees a
     * clean t_M (per-cell work) and t_N (per-particle insertion) split.
     * Spec §7: t_M = clear + prefix + offset-clear; t_N = cell_idx +
     * count+S + scatter. */
    struct timespec t0, t1;
    double t_M = 0, t_N = 0;
    int time_it = (t_M_out != NULL || t_N_out != NULL);

    /* 1. cell index per particle (clamped) -- O(N) */
    if (time_it) clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int i = 0; i < n; i++) {
        int cx = (int)floor(px[i] * inv);
        int cy = (int)floor(py[i] * inv);
        if (cx < 0)   cx = 0;
        if (cx >= gw) cx = gw - 1;
        if (cy < 0)   cy = 0;
        if (cy >= gh) cy = gh - 1;
        g->cell_idx[i] = cy * gw + cx;
    }
    if (time_it) { clock_gettime(CLOCK_MONOTONIC, &t1); t_N += ts_diff_s(&t0, &t1); }

    /* 2a. clear -- O(M) */
    if (time_it) clock_gettime(CLOCK_MONOTONIC, &t0);
    memset(g->cell_count, 0, sizeof(int) * (size_t)g->n_cells);
    if (time_it) { clock_gettime(CLOCK_MONOTONIC, &t1); t_M += ts_diff_s(&t0, &t1); }

    /* 2b. count + incremental S = sum n_i^2 -- O(N).
     *     Per particle: S += 2*count[cell]+1 BEFORE incrementing count.
     *     Adding one ball to a bin of size n raises n^2 by 2n+1.
     *     Exact S in O(N); no separate O(M) reduction over cell_count.    */
    if (time_it) clock_gettime(CLOCK_MONOTONIC, &t0);
    long S_acc = 0;
    for (int i = 0; i < n; i++) {
        int c = g->cell_idx[i];
        S_acc += 2L * (long)g->cell_count[c] + 1L;
        g->cell_count[c]++;
    }
    g->S_cached = S_acc;
    if (time_it) { clock_gettime(CLOCK_MONOTONIC, &t1); t_N += ts_diff_s(&t0, &t1); }

    /* 3. exclusive prefix sum -- O(M) */
    if (time_it) clock_gettime(CLOCK_MONOTONIC, &t0);
    int run = 0;
    for (int c = 0; c < g->n_cells; c++) {
        g->cell_start[c] = run;
        run += g->cell_count[c];
    }
    g->cell_start[g->n_cells] = run;
    if (time_it) { clock_gettime(CLOCK_MONOTONIC, &t1); t_M += ts_diff_s(&t0, &t1); }

    /* 4a. offset clear -- O(M) */
    int *offset = g->cell_count;
    if (time_it) clock_gettime(CLOCK_MONOTONIC, &t0);
    memset(offset, 0, sizeof(int) * (size_t)g->n_cells);
    if (time_it) { clock_gettime(CLOCK_MONOTONIC, &t1); t_M += ts_diff_s(&t0, &t1); }

    /* 4b. scatter -- O(N). After the loop offset[c] equals the original
     *     cell_count[c], so cell_count is naturally restored.            */
    if (time_it) clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int i = 0; i < n; i++) {
        int c = g->cell_idx[i];
        g->sorted_idx[g->cell_start[c] + offset[c]] = i;
        offset[c]++;
    }
    if (time_it) { clock_gettime(CLOCK_MONOTONIC, &t1); t_N += ts_diff_s(&t0, &t1); }

    if (t_M_out) *t_M_out = t_M;
    if (t_N_out) *t_N_out = t_N;
}

int  grid_M(const grid_t *g) { return g->n_cells; }

long grid_S(const grid_t *g) { return g->S_cached; }

long grid_count_S_at(double domain_w, double domain_h, int n,
                     const double *px, const double *py, double ell) {
    if (ell <= 0.0) return 0;
    int gw = (int)ceil(domain_w / ell);
    int gh = (int)ceil(domain_h / ell);
    if (gw < 1) gw = 1;
    if (gh < 1) gh = 1;
    int n_cells = gw * gh;
    int *counts = (int *)calloc((size_t)n_cells, sizeof(int));
    if (!counts) return -1;

    double inv = 1.0 / ell;
    long S_acc = 0;
    for (int i = 0; i < n; i++) {
        int cx = (int)floor(px[i] * inv);
        int cy = (int)floor(py[i] * inv);
        if (cx < 0)   cx = 0;
        if (cx >= gw) cx = gw - 1;
        if (cy < 0)   cy = 0;
        if (cy >= gh) cy = gh - 1;
        int c = cy * gw + cx;
        S_acc += 2L * (long)counts[c] + 1L;
        counts[c]++;
    }
    free(counts);
    return S_acc;
}

/* emit one candidate pair (ia, ib) in normalised order (i<j). */
static inline int emit_pair(int ia, int ib,
                             pair_t *out, int max_pairs, int *count) {
    if (*count >= max_pairs) return -1;
    int i = ia < ib ? ia : ib;
    int j = ia < ib ? ib : ia;
    out[*count].i = i;
    out[*count].j = j;
    (*count)++;
    return 0;
}

int grid_broad_phase(const grid_t *g,
                      pair_t *pairs_out, int max_pairs, int *count_out) {
    *count_out = 0;
    const int gw = g->gw, gh = g->gh;

    /* Half-stencil: (+1,0), (-1,+1), (0,+1), (+1,+1).
     * Each unordered pair of adjacent cells is visited from exactly one side. */
    static const int dxs[4] = {  1, -1,  0,  1 };
    static const int dys[4] = {  0,  1,  1,  1 };

    for (int cy = 0; cy < gh; cy++) {
        for (int cx = 0; cx < gw; cx++) {
            int c  = cy * gw + cx;
            int a0 = g->cell_start[c];
            int a1 = g->cell_start[c + 1];

            /* intra-cell pairs */
            for (int a = a0; a < a1; a++) {
                int ia = g->sorted_idx[a];
                for (int b = a + 1; b < a1; b++) {
                    int ib = g->sorted_idx[b];
                    if (emit_pair(ia, ib, pairs_out, max_pairs, count_out) != 0)
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
                        if (emit_pair(ia, ib, pairs_out, max_pairs, count_out) != 0)
                            return -1;
                    }
                }
            }
        }
    }
    return 0;
}

void narrow_phase(const pair_t *cands, int n_cands,
                   const double *px, const double *py, double r,
                   pair_t *out, int *n_overlap) {
    const double rr = (2.0 * r) * (2.0 * r);
    int w = 0;
    /* Safe with out == cands: writes always trail reads (w <= k). */
    for (int k = 0; k < n_cands; k++) {
        int i = cands[k].i, j = cands[k].j;
        double dx = px[j] - px[i];
        double dy = py[j] - py[i];
        if (dx * dx + dy * dy < rr) {
            out[w++] = cands[k];
        }
    }
    *n_overlap = w;
}

int grid_find_overlapping_pairs(const grid_t *g,
                                 const double *px, const double *py,
                                 double r,
                                 pair_t *pairs_out, int max_pairs, int *count_out) {
    int rc = grid_broad_phase(g, pairs_out, max_pairs, count_out);
    if (rc != 0) return rc;
    narrow_phase(pairs_out, *count_out, px, py, r, pairs_out, count_out);
    return 0;
}

double grid_broad_phase_iter_time(const grid_t *g, long *checksum_out) {
    /* Pure outer-iteration cost calibration: same cell-walk + half-stencil
     * neighbor structure as grid_broad_phase but no inner pair loops. The
     * actual M-prop cost is somewhere BETWEEN this (cell-walk only) and a
     * full broad-phase dryrun (which would also pick up the per-pair
     * for-loop overhead). We take the cell-walk-only flavour because
     *   (i)  per-pair overhead genuinely scales with S, so it belongs in
     *        t_S along with the actual emit work;
     *   (ii) this dryrun is cheap (one walk over cell_start), keeping the
     *        controller's overhead claim honest.
     * The consequence is a small residual M-prop bias in b_hat at low
     * occupancy, which makes the closed-form formula slightly undershoot
     * the per-frame oracle's ell. Documented in the Phase 5 summary. */
    const int gw = g->gw, gh = g->gh;
    static const int dxs[4] = {  1, -1,  0,  1 };
    static const int dys[4] = {  0,  1,  1,  1 };
    struct timespec ta, tb;
    long dummy = 0;
    clock_gettime(CLOCK_MONOTONIC, &ta);
    for (int cy = 0; cy < gh; cy++) {
        for (int cx = 0; cx < gw; cx++) {
            int c  = cy * gw + cx;
            int a0 = g->cell_start[c];
            int a1 = g->cell_start[c + 1];
            dummy += (long)(a1 - a0);
            for (int k = 0; k < 4; k++) {
                int nx = cx + dxs[k];
                int ny = cy + dys[k];
                if (nx < 0 || nx >= gw || ny < 0 || ny >= gh) continue;
                int nc = ny * gw + nx;
                int b0 = g->cell_start[nc];
                int b1 = g->cell_start[nc + 1];
                dummy += (long)(b1 - b0);
            }
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &tb);
    if (checksum_out) *checksum_out = dummy;
    return ts_diff_s(&ta, &tb);
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
