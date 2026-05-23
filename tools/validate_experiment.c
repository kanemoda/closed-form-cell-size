/*
 * Experiment-scale gate re-validation (spec §9.6 / TODO.md item 1).
 *
 * The Phase 1–3 unit gates run on a small std_cfg (N=400, 80×80). This tool
 * re-checks their essential PROPERTIES at the actual campaign scale of each
 * configs/experiment/*.json, before any campaign compute is spent:
 *
 *   - correctness : the grid broad+narrow overlap set == O(N²) brute force on a
 *                   mid-run frame (catches indexing/overflow at large M/N), 2D & 3D;
 *   - physics     : KE finite (no NaN/Inf) over the whole run; for force-free
 *                   scenarios KE is non-increasing (e<=1 dissipates);
 *   - determinism : two independent runs end bit-identical.
 *
 * Usage: ./validate_experiment configs/experiment/*.json
 */

#define _POSIX_C_SOURCE 199309L
#include "config.h"
#include "sim.h"
#include "scenarios.h"
#include "grid.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int cmp_pair(const void *a, const void *b) {
    const pair_t *pa = a, *pb = b;
    if (pa->i != pb->i) return pa->i - pb->i;
    return pa->j - pb->j;
}
static int pairsets_equal(pair_t *A, int na, pair_t *B, int nb) {
    if (na != nb) return 0;
    qsort(A, na, sizeof(pair_t), cmp_pair);
    qsort(B, nb, sizeof(pair_t), cmp_pair);
    for (int k = 0; k < na; k++)
        if (A[k].i != B[k].i || A[k].j != B[k].j) return 0;
    return 1;
}

/* Run the scenario to its final frame; return final KE and a position checksum
 * (for the determinism check). Tracks max KE and any non-finite value. */
static double run_to_end(const config_t *cfg, double *max_ke, int *finite_ok,
                         int *nincr_violations, double *checksum) {
    sim_t s; sim_init(&s, cfg);
    scenario_init(&s, scenario_from_name(cfg->scenario));
    int has_forces = scenario_has_forces(scenario_from_name(cfg->scenario));
    double ke_prev = sim_total_kinetic_energy(&s);
    *max_ke = ke_prev; *finite_ok = 1; *nincr_violations = 0;
    for (int t = 0; t < cfg->num_frames; t++) {
        sim_step(&s, NULL);
        double ke = sim_total_kinetic_energy(&s);
        if (!isfinite(ke)) *finite_ok = 0;
        if (ke > *max_ke) *max_ke = ke;
        if (!has_forces && ke > ke_prev * (1.0 + 1e-9) + 1e-12) (*nincr_violations)++;
        ke_prev = ke;
    }
    double cs = 0.0;
    for (int i = 0; i < s.n; i++) {
        cs += s.px[i] * 1.000001 + s.py[i] * 1.0000003;
        if (cfg->dim == 3) cs += s.pz[i] * 0.9999997;
    }
    double ke_final = ke_prev;
    *checksum = cs;
    sim_free(&s);
    return ke_final;
}

static int validate_one(const char *path) {
    config_t cfg;
    if (config_load(path, &cfg) != 0) { fprintf(stderr, "  %-28s LOAD FAILED\n", path); return 0; }
    cfg.log_enabled = 0; cfg.log_csv_path[0] = '\0';
    cfg.snapshot_csv_path[0] = '\0'; cfg.oracle_csv_path[0] = '\0'; cfg.adapter_mode[0] = '\0';
    int has_forces = scenario_has_forces(scenario_from_name(cfg.scenario));

    /* --- correctness: grid == brute on a mid-run frame --- */
    sim_t s; sim_init(&s, &cfg);
    scenario_init(&s, scenario_from_name(cfg.scenario));
    int probe_frame = cfg.num_frames / 2;
    for (int t = 0; t <= probe_frame; t++) sim_step(&s, NULL);
    long cap = 64L * cfg.N;
    pair_t *gp = malloc(sizeof(pair_t) * (size_t)cap);
    pair_t *bp = malloc(sizeof(pair_t) * (size_t)cap);
    grid_t g;
    int ng = 0, nb = 0, correctness_ok = 1;
    if (cfg.dim == 3) {
        grid_init3d(&g, s.n, cfg.domain_w, cfg.domain_h, cfg.domain_d, cfg.cell_size);
        grid_build3d(&g, s.px, s.py, s.pz);
        if (grid_find_overlapping_pairs3d(&g, s.px, s.py, s.pz, cfg.radius, gp, (int)cap, &ng) != 0) correctness_ok = 0;
        bruteforce_find_overlapping_pairs3d(s.n, s.px, s.py, s.pz, cfg.radius, bp, (int)cap, &nb);
    } else {
        grid_init(&g, s.n, cfg.domain_w, cfg.domain_h, cfg.cell_size);
        grid_build(&g, s.px, s.py);
        if (grid_find_overlapping_pairs(&g, s.px, s.py, cfg.radius, gp, (int)cap, &ng) != 0) correctness_ok = 0;
        bruteforce_find_overlapping_pairs(s.n, s.px, s.py, cfg.radius, bp, (int)cap, &nb);
    }
    if (!pairsets_equal(gp, ng, bp, nb)) correctness_ok = 0;
    grid_free(&g); free(gp); free(bp);
    sim_free(&s);

    /* --- physics + determinism: two full runs --- */
    double max_ke1, max_ke2, cs1, cs2, kef1, kef2;
    int fin1, fin2, nincr1, nincr2;
    kef1 = run_to_end(&cfg, &max_ke1, &fin1, &nincr1, &cs1);
    kef2 = run_to_end(&cfg, &max_ke2, &fin2, &nincr2, &cs2);

    int finite_ok   = fin1 && fin2;
    int det_ok      = (kef1 == kef2) && (cs1 == cs2);
    int nincr_ok    = has_forces ? 1 : (nincr1 == 0);
    int ok = correctness_ok && finite_ok && det_ok && nincr_ok;

    const char *base = strrchr(path, '/'); base = base ? base + 1 : path;
    fprintf(stderr,
        "  %-22s N=%-6d %s | overlaps=%-6d grid==brute:%s | KE0->max=%.3g->%.3g fin:%s %s:%s | det:%s | %s\n",
        base, cfg.N, cfg.dim == 3 ? "3D" : "2D",
        nb, correctness_ok ? "Y" : "N",
        kef1, max_ke1, finite_ok ? "Y" : "N",
        has_forces ? "force" : "KE↓", has_forces ? "--" : (nincr_ok ? "Y" : "N"),
        det_ok ? "Y" : "N", ok ? "OK" : "FAIL");
    return ok;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s config.json [config.json ...]\n", argv[0]); return 2; }
    fprintf(stderr, "=== Experiment-scale gate re-validation (spec §9.6) ===\n");
    int all_ok = 1, n = 0;
    for (int i = 1; i < argc; i++) { all_ok &= validate_one(argv[i]); n++; }
    fprintf(stderr, "\n%s  (%d configs)\n", all_ok ? "ALL OK" : "SOME FAILED", n);
    return all_ok ? 0 : 1;
}
