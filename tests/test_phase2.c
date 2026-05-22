/*
 * Phase 2 verification gates (spec section 10):
 *   G1. S logged == independent recompute of sum(n_i^2), every frame.
 *   G2. num_cells == grid.n_cells, and == round(V / ell^d) when W,H are
 *       integer multiples of ell.
 *   G3. broad+narrow pipeline output == brute force overlap set;
 *       candidate_pairs is a superset of collisions.
 *   G4. All four timers are positive; t_grid + t_broad + t_narrow <= t_total.
 *   G5. Two runs of the same config produce a CSV whose deterministic
 *       columns are bit-identical (timer columns are wall-clock and not
 *       expected to match -- they are checked separately).
 */

#include "rng.h"
#include "config.h"
#include "sim.h"
#include "grid.h"
#include "log.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---------------------------------------------------------------- helpers */

/* In-place CSV tokenizer. Returns the next comma- or newline-delimited
 * field, advances *s past it, and writes '\0' at the delimiter. */
static char *csv_next(char **s) {
    if (!*s || !**s) return NULL;
    char *start = *s;
    while (**s && **s != ',' && **s != '\n') (*s)++;
    if (**s) { **s = '\0'; (*s)++; }
    return start;
}

static int cmp_pair(const void *a, const void *b) {
    const pair_t *pa = a, *pb = b;
    if (pa->i != pb->i) return pa->i - pb->i;
    return pa->j - pb->j;
}

/* ----------------------------------------- G1: S recompute == logged S */
static int g1_S_recompute(void) {
    config_t cfg; config_default(&cfg);
    cfg.N = 500; cfg.domain_w = 80; cfg.domain_h = 80;
    cfg.radius = 0.5; cfg.cell_size = 2.0; cfg.dt = 0.01;
    cfg.init_speed = 5; cfg.num_frames = 50; cfg.seed = 1;
    cfg.restitution = 0.9;

    sim_t s; sim_init(&s, &cfg); sim_init_random(&s);

    int ok = 1, bad = 0;
    long worst_diff = 0;
    for (int t = 0; t < cfg.num_frames; t++) {
        sim_metrics_t m = { .frame = t };
        sim_step(&s, &m);

        /* Independently recompute S from grid contents.
         * (Same formula, but a separate, post-step recount.)            */
        long S_recompute = 0;
        for (int c = 0; c < s.grid.n_cells; c++) {
            long ni = s.grid.cell_count[c];
            S_recompute += ni * ni;
        }
        if (S_recompute != m.S) {
            long d = S_recompute - m.S; if (d < 0) d = -d;
            if (d > worst_diff) worst_diff = d;
            bad++;
            ok = 0;
        }
    }
    fprintf(stderr, "  G1 S recompute over %d frames: bad=%d worst_diff=%ld  %s\n",
            cfg.num_frames, bad, worst_diff, ok ? "OK" : "FAIL");
    sim_free(&s);
    return ok;
}

/* ------------------------------------ G2: num_cells == grid && V/ell^d */
struct g2_case { int N; double W, H, cell; };
static int g2_num_cells(void) {
    /* All cases here have W and H as integer multiples of cell, so the
     * ceil-based count equals round(V / ell^d). */
    struct g2_case cases[] = {
        { 200,  80,  80, 2.0 },     /* 40x40 = 1600 cells; V/ell^d = 1600 */
        { 200, 100,  60, 2.0 },     /* 50x30 = 1500;       V/ell^d = 1500 */
        { 500, 120, 120, 1.5 },     /* 80x80 = 6400;       V/ell^d = 6400 */
        { 300,  50,  50, 1.0 },     /* 50x50 = 2500;       V/ell^d = 2500 */
        { 100,  20,  10, 1.0 },     /* 20x10 = 200;        V/ell^d = 200  */
    };
    const int n_cases = (int)(sizeof(cases) / sizeof(cases[0]));

    int ok = 1;
    for (int k = 0; k < n_cases; k++) {
        config_t cfg; config_default(&cfg);
        cfg.N = cases[k].N; cfg.domain_w = cases[k].W; cfg.domain_h = cases[k].H;
        cfg.radius = 0.4; cfg.cell_size = cases[k].cell; cfg.dt = 0.01;
        cfg.init_speed = 1; cfg.num_frames = 1; cfg.seed = (uint64_t)(k + 1);

        sim_t s; sim_init(&s, &cfg); sim_init_random(&s);
        sim_metrics_t m = { .frame = 0 };
        sim_step(&s, &m);

        int expected_grid  = (int)ceil(cfg.domain_w / cfg.cell_size)
                           * (int)ceil(cfg.domain_h / cfg.cell_size);
        int expected_model = (int)lround((cfg.domain_w * cfg.domain_h) /
                                         (cfg.cell_size * cfg.cell_size));
        int case_ok = (m.num_cells == s.grid.n_cells)
                   && (m.num_cells == expected_grid)
                   && (m.num_cells == expected_model);
        if (!case_ok) ok = 0;
        fprintf(stderr,
            "  G2 case W=%g H=%g cell=%g: num_cells=%d grid=%d V/ell^d=%d  %s\n",
            cfg.domain_w, cfg.domain_h, cfg.cell_size,
            m.num_cells, expected_grid, expected_model,
            case_ok ? "ok" : "MISMATCH");
        sim_free(&s);
    }
    fprintf(stderr, "  G2 num_cells: %s\n", ok ? "OK" : "FAIL");
    return ok;
}

/* ----------- G3: pipeline output == brute force; candidates superset */
static int g3_pipeline_eq_brute(void) {
    /* Dense, overlapping scene -- random placement with no rejection. */
    int    N    = 400;
    double W    = 15, H = 15;
    double r    = 0.5, cell = 1.0;
    uint64_t sd = 31;

    rng_state_t rng; rng_seed(&rng, sd);
    double *px = (double*)malloc(sizeof(double)*N);
    double *py = (double*)malloc(sizeof(double)*N);
    for (int i = 0; i < N; i++) {
        px[i] = rng_uniform(&rng, r, W - r);
        py[i] = rng_uniform(&rng, r, H - r);
    }

    grid_t g; grid_init(&g, N, W, H, cell); grid_build(&g, px, py);

    int max_pairs = N * (N - 1) / 2 + 1;
    pair_t *cands  = (pair_t*)malloc(sizeof(pair_t) * (size_t)max_pairs);
    pair_t *pipe   = (pair_t*)malloc(sizeof(pair_t) * (size_t)max_pairs);
    pair_t *brute  = (pair_t*)malloc(sizeof(pair_t) * (size_t)max_pairs);

    int n_cand = 0, n_pipe = 0, n_brute = 0;
    int rc1 = grid_broad_phase(&g, cands, max_pairs, &n_cand);
    narrow_phase(cands, n_cand, px, py, r, pipe, &n_pipe);
    int rc2 = bruteforce_find_overlapping_pairs(N, px, py, r, brute, max_pairs, &n_brute);

    qsort(pipe,  n_pipe,  sizeof(pair_t), cmp_pair);
    qsort(brute, n_brute, sizeof(pair_t), cmp_pair);

    int set_eq = (n_pipe == n_brute);
    if (set_eq) {
        for (int k = 0; k < n_pipe; k++) {
            if (pipe[k].i != brute[k].i || pipe[k].j != brute[k].j) {
                set_eq = 0; break;
            }
        }
    }

    /* Superset check: every overlap pair must appear in the candidate list. */
    qsort(cands, n_cand, sizeof(pair_t), cmp_pair);
    int superset = 1, ci = 0;
    for (int k = 0; k < n_brute; k++) {
        while (ci < n_cand &&
               (cands[ci].i < brute[k].i ||
                (cands[ci].i == brute[k].i && cands[ci].j < brute[k].j))) ci++;
        if (ci >= n_cand || cands[ci].i != brute[k].i || cands[ci].j != brute[k].j) {
            superset = 0; break;
        }
    }

    int count_invariant = n_cand >= n_pipe;

    int ok = (rc1 == 0) && (rc2 == 0) && set_eq && superset && count_invariant;
    fprintf(stderr,
        "  G3 N=%d W=%g cell=%g: cands=%d pipe=%d brute=%d  "
        "pipe==brute=%d cands_superset=%d cands>=pipe=%d  %s\n",
        N, W, cell, n_cand, n_pipe, n_brute,
        set_eq, superset, count_invariant, ok ? "OK" : "FAIL");

    free(cands); free(pipe); free(brute);
    free(px); free(py);
    grid_free(&g);
    return ok;
}

/* ----- G4: positive timers and t_grid + t_broad + t_narrow <= t_total */
static int g4_timer_invariants(void) {
    config_t cfg; config_default(&cfg);
    cfg.N = 1000; cfg.domain_w = 100; cfg.domain_h = 100; cfg.radius = 0.5;
    cfg.cell_size = 2.0; cfg.dt = 0.01; cfg.init_speed = 5;
    cfg.num_frames = 30; cfg.seed = 11; cfg.restitution = 0.9;

    sim_t s; sim_init(&s, &cfg); sim_init_random(&s);

    int ok = 1;
    int bad_pos = 0, bad_inv = 0;
    double worst_excess = 0.0;
    for (int t = 0; t < cfg.num_frames; t++) {
        sim_metrics_t m = { .frame = t };
        sim_step(&s, &m);
        if (m.t_grid   <= 0.0 ||
            m.t_broad  <= 0.0 ||
            m.t_narrow <= 0.0 ||
            m.t_total  <= 0.0) { bad_pos++; ok = 0; }
        double sub = m.t_grid + m.t_broad + m.t_narrow;
        if (sub > m.t_total) {
            double excess = sub - m.t_total;
            if (excess > worst_excess) worst_excess = excess;
            bad_inv++;
            ok = 0;
        }
    }
    fprintf(stderr,
        "  G4 timer invariants (%d frames): bad_positive=%d bad_inequality=%d worst_excess=%g  %s\n",
        cfg.num_frames, bad_pos, bad_inv, worst_excess, ok ? "OK" : "FAIL");
    sim_free(&s);
    return ok;
}

/* ---------- G5: deterministic CSV columns bit-identical across runs */
static int g5_csv_determinism(void) {
    const char *path1 = "/tmp/phase2_det_run1.csv";
    const char *path2 = "/tmp/phase2_det_run2.csv";

    config_t cfg; config_default(&cfg);
    cfg.N = 200; cfg.domain_w = 30; cfg.domain_h = 30; cfg.radius = 0.5;
    cfg.cell_size = 1.5; cfg.dt = 0.01; cfg.init_speed = 4;
    cfg.num_frames = 40; cfg.seed = 12345; cfg.restitution = 0.9;

    for (int run = 0; run < 2; run++) {
        sim_t s; sim_init(&s, &cfg); sim_init_random(&s);
        FILE *csv = fopen(run == 0 ? path1 : path2, "w");
        log_write_header(csv);
        for (int t = 0; t < cfg.num_frames; t++) {
            sim_metrics_t m = { .frame = t };
            sim_step(&s, &m);
            log_write_row(csv, &m);
        }
        fclose(csv);
        sim_free(&s);
    }

    FILE *f1 = fopen(path1, "r"), *f2 = fopen(path2, "r");
    int ok = 1, row = 0, det_mismatch = 0, timer_violations = 0;
    char line1[1024], line2[1024];

    /* Header line must match exactly. */
    if (!fgets(line1, sizeof(line1), f1) || !fgets(line2, sizeof(line2), f2) ||
        strcmp(line1, line2) != 0) {
        fprintf(stderr, "  G5 header mismatch\n");
        ok = 0;
    }
    row++;

    /* Column positions:
     *   0 frame, 1 cell_size, 2 N, 3 num_cells, 4 S,
     *   5 candidate_pairs, 6 collisions,
     *   7 t_grid, 8 t_broad, 9 t_narrow, 10 t_total,
     *   11 kinetic_energy, 12 momentum_magnitude
     */
    const int det_cols[]   = { 0, 1, 2, 3, 4, 5, 6, 11, 12 };
    const int n_det        = (int)(sizeof(det_cols) / sizeof(det_cols[0]));
    const int timer_cols[] = { 7, 8, 9, 10 };

    while (fgets(line1, sizeof(line1), f1) && fgets(line2, sizeof(line2), f2)) {
        char buf1[1024], buf2[1024];
        memcpy(buf1, line1, sizeof(buf1));
        memcpy(buf2, line2, sizeof(buf2));
        char *s1 = buf1, *s2 = buf2;
        char *toks1[13] = {0}, *toks2[13] = {0};
        for (int c = 0; c < 13; c++) {
            toks1[c] = csv_next(&s1);
            toks2[c] = csv_next(&s2);
            if (!toks1[c] || !toks2[c]) { ok = 0; goto done; }
        }
        for (int k = 0; k < n_det; k++) {
            int c = det_cols[k];
            if (strcmp(toks1[c], toks2[c]) != 0) {
                if (det_mismatch < 3) {
                    fprintf(stderr,
                        "  G5 row %d col %d differs: '%s' vs '%s'\n",
                        row, c, toks1[c], toks2[c]);
                }
                det_mismatch++;
                ok = 0;
            }
        }
        /* Timer columns: just check positive and sum inequality, per run 1. */
        double tg = atof(toks1[timer_cols[0]]);
        double tb = atof(toks1[timer_cols[1]]);
        double tn = atof(toks1[timer_cols[2]]);
        double tt = atof(toks1[timer_cols[3]]);
        if (tg <= 0 || tb <= 0 || tn <= 0 || tt <= 0) timer_violations++;
        if (tg + tb + tn > tt)                       timer_violations++;
        row++;
    }
    if (timer_violations) { ok = 0; }
done:
    fclose(f1); fclose(f2);

    fprintf(stderr,
        "  G5 CSV determinism (%d data rows): det_mismatches=%d  timer_violations=%d  %s\n",
        row - 1, det_mismatch, timer_violations, ok ? "OK" : "FAIL");
    return ok;
}

int main(void) {
    fprintf(stderr, "=== Phase 2 Verification Gates ===\n");
    int all_ok = 1;
    all_ok &= g1_S_recompute();
    all_ok &= g2_num_cells();
    all_ok &= g3_pipeline_eq_brute();
    all_ok &= g4_timer_invariants();
    all_ok &= g5_csv_determinism();
    fprintf(stderr, "\n%s\n", all_ok ? "Phase 2: ALL OK" : "Phase 2: FAILED");
    return all_ok ? 0 : 1;
}
