#define _POSIX_C_SOURCE 199309L
#include "oracle.h"
#include "grid.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static inline double ts_diff_s(const struct timespec *a, const struct timespec *b) {
    return (double)(b->tv_sec - a->tv_sec) + (double)(b->tv_nsec - a->tv_nsec) * 1e-9;
}

static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    if (da < db) return -1;
    if (da > db) return  1;
    return 0;
}

static double median_double(double *v, int n) {
    qsort(v, n, sizeof(double), cmp_double);
    if (n & 1) return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

int oracle_make_geo_sweep(double ell_min, double ell_max, double gamma,
                          double *ells, int max_pts) {
    if (gamma <= 1.0 || ell_min <= 0.0 || ell_max < ell_min) return 0;
    int n = 0;
    double cur = ell_min;
    while (cur <= ell_max && n < max_pts) {
        ells[n++] = cur;
        cur *= gamma;
    }
    return n;
}

void oracle_eval(int n, const double *px, const double *py, double radius,
                 double domain_w, double domain_h,
                 double ell, int K,
                 oracle_point_t *out) {
    grid_t g;
    grid_init(&g, n, domain_w, domain_h, ell);

    int pairs_cap = 16 * n;
    if (pairs_cap < 1024) pairs_cap = 1024;
    pair_t *pairs = (pair_t *)malloc(sizeof(pair_t) * (size_t)pairs_cap);
    if (!pairs) { fprintf(stderr, "oracle_eval: malloc failed\n"); exit(1); }

    /* K must be >= 1; clamp. */
    if (K < 1) K = 1;
    double *t_g  = (double *)malloc(sizeof(double) * (size_t)K);
    double *t_b  = (double *)malloc(sizeof(double) * (size_t)K);
    double *t_nw = (double *)malloc(sizeof(double) * (size_t)K);

    int n_cand = 0;
    long S_last = 0;
    int  M_last = 0;

    for (int k = 0; k < K; k++) {
        struct timespec ta, tb, tc, td;
        clock_gettime(CLOCK_MONOTONIC, &ta);
        grid_build(&g, px, py);
        clock_gettime(CLOCK_MONOTONIC, &tb);

        while (grid_broad_phase(&g, pairs, pairs_cap, &n_cand) != 0) {
            pairs_cap *= 2;
            pairs = (pair_t *)realloc(pairs, sizeof(pair_t) * (size_t)pairs_cap);
            if (!pairs) { fprintf(stderr, "oracle_eval: realloc failed\n"); exit(1); }
        }
        clock_gettime(CLOCK_MONOTONIC, &tc);

        int n_overlap;
        narrow_phase(pairs, n_cand, px, py, radius, pairs, &n_overlap);
        clock_gettime(CLOCK_MONOTONIC, &td);

        t_g [k] = ts_diff_s(&ta, &tb);
        t_b [k] = ts_diff_s(&tb, &tc);
        t_nw[k] = ts_diff_s(&tc, &td);

        S_last = grid_S(&g);
        M_last = grid_M(&g);
    }

    out->ell             = ell;
    out->num_cells       = M_last;
    out->S               = S_last;
    out->candidate_pairs = n_cand;
    out->t_grid          = median_double(t_g,  K);
    out->t_broad         = median_double(t_b,  K);
    out->t_narrow        = median_double(t_nw, K);
    out->t_detect        = out->t_grid + out->t_broad + out->t_narrow;

    free(t_g); free(t_b); free(t_nw);
    free(pairs);
    grid_free(&g);
}

void oracle_eval_sweep(int n, const double *px, const double *py, double radius,
                       double domain_w, double domain_h,
                       const double *ells, int n_ells, int K,
                       oracle_point_t *out) {
    for (int i = 0; i < n_ells; i++) {
        oracle_eval(n, px, py, radius, domain_w, domain_h,
                    ells[i], K, &out[i]);
    }
}

void oracle_eval_split(int n, const double *px, const double *py, double radius,
                       double domain_w, double domain_h,
                       double ell, int K,
                       oracle_point_t *out_sel, oracle_point_t *out_eval) {
    grid_t g;
    grid_init(&g, n, domain_w, domain_h, ell);

    int pairs_cap = 16 * n;
    if (pairs_cap < 1024) pairs_cap = 1024;
    pair_t *pairs = (pair_t *)malloc(sizeof(pair_t) * (size_t)pairs_cap);
    if (!pairs) { fprintf(stderr, "oracle_eval_split: malloc failed\n"); exit(1); }

    if (K < 1) K = 1;
    /* 2K full-pipeline timings, interleaved into SEL (even k) / EVAL (odd k)
     * so both groups see the same systematic conditions (no A-then-B drift). */
    double *g_sel = (double *)malloc(sizeof(double) * (size_t)K);
    double *b_sel = (double *)malloc(sizeof(double) * (size_t)K);
    double *n_sel = (double *)malloc(sizeof(double) * (size_t)K);
    double *g_ev  = (double *)malloc(sizeof(double) * (size_t)K);
    double *b_ev  = (double *)malloc(sizeof(double) * (size_t)K);
    double *n_ev  = (double *)malloc(sizeof(double) * (size_t)K);

    int n_cand = 0;
    long S_last = 0;
    int  M_last = 0;

    for (int k = 0; k < 2 * K; k++) {
        struct timespec ta, tb, tc, td;
        clock_gettime(CLOCK_MONOTONIC, &ta);
        grid_build(&g, px, py);
        clock_gettime(CLOCK_MONOTONIC, &tb);

        while (grid_broad_phase(&g, pairs, pairs_cap, &n_cand) != 0) {
            pairs_cap *= 2;
            pairs = (pair_t *)realloc(pairs, sizeof(pair_t) * (size_t)pairs_cap);
            if (!pairs) { fprintf(stderr, "oracle_eval_split: realloc failed\n"); exit(1); }
        }
        clock_gettime(CLOCK_MONOTONIC, &tc);

        int n_overlap;
        narrow_phase(pairs, n_cand, px, py, radius, pairs, &n_overlap);
        clock_gettime(CLOCK_MONOTONIC, &td);

        int idx = k / 2;
        if ((k & 1) == 0) {
            g_sel[idx] = ts_diff_s(&ta, &tb);
            b_sel[idx] = ts_diff_s(&tb, &tc);
            n_sel[idx] = ts_diff_s(&tc, &td);
        } else {
            g_ev[idx]  = ts_diff_s(&ta, &tb);
            b_ev[idx]  = ts_diff_s(&tb, &tc);
            n_ev[idx]  = ts_diff_s(&tc, &td);
        }
        S_last = grid_S(&g);
        M_last = grid_M(&g);
    }

    out_sel->ell  = ell;  out_eval->ell  = ell;
    out_sel->num_cells = M_last;  out_eval->num_cells = M_last;
    out_sel->S    = S_last;       out_eval->S    = S_last;
    out_sel->candidate_pairs = n_cand;  out_eval->candidate_pairs = n_cand;
    out_sel->t_grid   = median_double(g_sel, K);
    out_sel->t_broad  = median_double(b_sel, K);
    out_sel->t_narrow = median_double(n_sel, K);
    out_sel->t_detect = out_sel->t_grid + out_sel->t_broad + out_sel->t_narrow;
    out_eval->t_grid   = median_double(g_ev, K);
    out_eval->t_broad  = median_double(b_ev, K);
    out_eval->t_narrow = median_double(n_ev, K);
    out_eval->t_detect = out_eval->t_grid + out_eval->t_broad + out_eval->t_narrow;

    free(g_sel); free(b_sel); free(n_sel);
    free(g_ev);  free(b_ev);  free(n_ev);
    free(pairs);
    grid_free(&g);
}

/* ===================================================================== 3D */

void oracle_eval3d(int n, const double *px, const double *py, const double *pz,
                   double radius, double domain_w, double domain_h, double domain_d,
                   double ell, int K, oracle_point_t *out) {
    grid_t g;
    grid_init3d(&g, n, domain_w, domain_h, domain_d, ell);

    int pairs_cap = 16 * n;
    if (pairs_cap < 1024) pairs_cap = 1024;
    pair_t *pairs = (pair_t *)malloc(sizeof(pair_t) * (size_t)pairs_cap);
    if (!pairs) { fprintf(stderr, "oracle_eval3d: malloc failed\n"); exit(1); }

    if (K < 1) K = 1;
    double *t_g  = (double *)malloc(sizeof(double) * (size_t)K);
    double *t_b  = (double *)malloc(sizeof(double) * (size_t)K);
    double *t_nw = (double *)malloc(sizeof(double) * (size_t)K);

    int n_cand = 0;
    long S_last = 0;
    int  M_last = 0;
    for (int k = 0; k < K; k++) {
        struct timespec ta, tb, tc, td;
        clock_gettime(CLOCK_MONOTONIC, &ta);
        grid_build3d(&g, px, py, pz);
        clock_gettime(CLOCK_MONOTONIC, &tb);
        while (grid_broad_phase(&g, pairs, pairs_cap, &n_cand) != 0) {
            pairs_cap *= 2;
            pairs = (pair_t *)realloc(pairs, sizeof(pair_t) * (size_t)pairs_cap);
            if (!pairs) { fprintf(stderr, "oracle_eval3d: realloc failed\n"); exit(1); }
        }
        clock_gettime(CLOCK_MONOTONIC, &tc);
        int n_overlap;
        narrow_phase3d(pairs, n_cand, px, py, pz, radius, pairs, &n_overlap);
        clock_gettime(CLOCK_MONOTONIC, &td);
        t_g[k] = ts_diff_s(&ta, &tb);
        t_b[k] = ts_diff_s(&tb, &tc);
        t_nw[k] = ts_diff_s(&tc, &td);
        S_last = grid_S(&g);
        M_last = grid_M(&g);
    }

    out->ell = ell;
    out->num_cells = M_last;
    out->S = S_last;
    out->candidate_pairs = n_cand;
    out->t_grid   = median_double(t_g,  K);
    out->t_broad  = median_double(t_b,  K);
    out->t_narrow = median_double(t_nw, K);
    out->t_detect = out->t_grid + out->t_broad + out->t_narrow;

    free(t_g); free(t_b); free(t_nw);
    free(pairs);
    grid_free(&g);
}

void oracle_eval_sweep3d(int n, const double *px, const double *py, const double *pz,
                         double radius, double domain_w, double domain_h, double domain_d,
                         const double *ells, int n_ells, int K, oracle_point_t *out) {
    for (int i = 0; i < n_ells; i++)
        oracle_eval3d(n, px, py, pz, radius, domain_w, domain_h, domain_d,
                      ells[i], K, &out[i]);
}

void oracle_eval_split3d(int n, const double *px, const double *py, const double *pz,
                         double radius, double domain_w, double domain_h, double domain_d,
                         double ell, int K,
                         oracle_point_t *out_sel, oracle_point_t *out_eval) {
    grid_t g;
    grid_init3d(&g, n, domain_w, domain_h, domain_d, ell);

    int pairs_cap = 16 * n;
    if (pairs_cap < 1024) pairs_cap = 1024;
    pair_t *pairs = (pair_t *)malloc(sizeof(pair_t) * (size_t)pairs_cap);
    if (!pairs) { fprintf(stderr, "oracle_eval_split3d: malloc failed\n"); exit(1); }

    if (K < 1) K = 1;
    double *g_sel = (double *)malloc(sizeof(double) * (size_t)K);
    double *b_sel = (double *)malloc(sizeof(double) * (size_t)K);
    double *n_sel = (double *)malloc(sizeof(double) * (size_t)K);
    double *g_ev  = (double *)malloc(sizeof(double) * (size_t)K);
    double *b_ev  = (double *)malloc(sizeof(double) * (size_t)K);
    double *n_ev  = (double *)malloc(sizeof(double) * (size_t)K);

    int n_cand = 0;
    long S_last = 0;
    int  M_last = 0;
    for (int k = 0; k < 2 * K; k++) {
        struct timespec ta, tb, tc, td;
        clock_gettime(CLOCK_MONOTONIC, &ta);
        grid_build3d(&g, px, py, pz);
        clock_gettime(CLOCK_MONOTONIC, &tb);
        while (grid_broad_phase(&g, pairs, pairs_cap, &n_cand) != 0) {
            pairs_cap *= 2;
            pairs = (pair_t *)realloc(pairs, sizeof(pair_t) * (size_t)pairs_cap);
            if (!pairs) { fprintf(stderr, "oracle_eval_split3d: realloc failed\n"); exit(1); }
        }
        clock_gettime(CLOCK_MONOTONIC, &tc);
        int n_overlap;
        narrow_phase3d(pairs, n_cand, px, py, pz, radius, pairs, &n_overlap);
        clock_gettime(CLOCK_MONOTONIC, &td);
        int idx = k / 2;
        if ((k & 1) == 0) {
            g_sel[idx] = ts_diff_s(&ta, &tb);
            b_sel[idx] = ts_diff_s(&tb, &tc);
            n_sel[idx] = ts_diff_s(&tc, &td);
        } else {
            g_ev[idx]  = ts_diff_s(&ta, &tb);
            b_ev[idx]  = ts_diff_s(&tb, &tc);
            n_ev[idx]  = ts_diff_s(&tc, &td);
        }
        S_last = grid_S(&g);
        M_last = grid_M(&g);
    }

    out_sel->ell = ell;  out_eval->ell = ell;
    out_sel->num_cells = M_last;  out_eval->num_cells = M_last;
    out_sel->S = S_last;          out_eval->S = S_last;
    out_sel->candidate_pairs = n_cand;  out_eval->candidate_pairs = n_cand;
    out_sel->t_grid   = median_double(g_sel, K);
    out_sel->t_broad  = median_double(b_sel, K);
    out_sel->t_narrow = median_double(n_sel, K);
    out_sel->t_detect = out_sel->t_grid + out_sel->t_broad + out_sel->t_narrow;
    out_eval->t_grid   = median_double(g_ev, K);
    out_eval->t_broad  = median_double(b_ev, K);
    out_eval->t_narrow = median_double(n_ev, K);
    out_eval->t_detect = out_eval->t_grid + out_eval->t_broad + out_eval->t_narrow;

    free(g_sel); free(b_sel); free(n_sel);
    free(g_ev);  free(b_ev);  free(n_ev);
    free(pairs);
    grid_free(&g);
}
