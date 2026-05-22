#include "baselines.h"
#include "sim.h"
#include "scenarios.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

double baseline_uniform_analytical(int N, double domain_w, double domain_h,
                                   double radius) {
    double V    = domain_w * domain_h;
    double ell  = sqrt(V / (double)N);
    double flr  = 2.0 * radius;
    return ell > flr ? ell : flr;
}

void run_sim_with_cell_size(const config_t *cfg, double cell_size,
                            sim_run_summary_t *out) {
    config_t mod = *cfg;
    mod.cell_size   = cell_size;
    mod.log_enabled = 0;                 /* no CSV during oracle runs */
    mod.log_csv_path[0]      = '\0';
    mod.snapshot_csv_path[0] = '\0';

    sim_t s; sim_init(&s, &mod);
    scenario_id_t id = scenario_from_name(mod.scenario);
    scenario_init(&s, id);

    double sum_g = 0, sum_b = 0, sum_nw = 0;
    long   tot_c = 0, tot_oc = 0;
    int    frames = mod.num_frames;
    for (int t = 0; t < frames; t++) {
        sim_metrics_t m = { .frame = t };
        sim_step(&s, &m);
        sum_g  += m.t_grid;
        sum_b  += m.t_broad;
        sum_nw += m.t_narrow;
        tot_c  += m.candidate_pairs;
        tot_oc += m.collisions;
    }
    out->mean_t_grid     = sum_g  / (double)frames;
    out->mean_t_broad    = sum_b  / (double)frames;
    out->mean_t_narrow   = sum_nw / (double)frames;
    out->mean_t_detect   = (sum_g + sum_b + sum_nw) / (double)frames;
    out->total_candidates = tot_c;
    out->total_collisions = tot_oc;
    out->frames           = frames;

    sim_free(&s);
}

void static_oracle_run(const config_t *cfg, const double *ells, int n_ells,
                       static_oracle_result_t *out) {
    out->n_candidates = n_ells;
    out->ells         = (double *)malloc(sizeof(double) * (size_t)n_ells);
    out->mean_detects = (double *)malloc(sizeof(double) * (size_t)n_ells);
    memcpy(out->ells, ells, sizeof(double) * (size_t)n_ells);

    out->best_index       = 0;
    out->best_mean_detect = 1.0e18;
    for (int i = 0; i < n_ells; i++) {
        sim_run_summary_t s;
        run_sim_with_cell_size(cfg, ells[i], &s);
        out->mean_detects[i] = s.mean_t_detect;
        if (s.mean_t_detect < out->best_mean_detect) {
            out->best_mean_detect = s.mean_t_detect;
            out->best_index       = i;
        }
    }
    out->best_ell = out->ells[out->best_index];
}

void static_oracle_free(static_oracle_result_t *r) {
    free(r->ells);
    free(r->mean_detects);
    r->ells         = NULL;
    r->mean_detects = NULL;
}
