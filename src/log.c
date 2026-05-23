#include "log.h"
#include <stdio.h>

FILE *log_open(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "log_open: failed to open '%s' for writing\n", path);
        return NULL;
    }
    return f;
}

void log_write_header(FILE *f) {
    if (!f) return;
    fputs("frame,cell_size,N,num_cells,S,candidate_pairs,collisions,"
          "t_grid,t_broad,t_narrow,t_total,t_M,t_N,"
          "kinetic_energy,momentum_magnitude\n", f);
}

void log_write_row(FILE *f, const sim_metrics_t *m) {
    if (!f || !m) return;
    fprintf(f,
            "%d,%.17g,%d,%d,%ld,%d,%d,"          /* deterministic head */
            "%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,"     /* wall-clock timers */
            "%.17g,%.17g\n",                     /* deterministic tail */
            m->frame, m->cell_size, m->N, m->num_cells, m->S,
            m->candidate_pairs, m->collisions,
            m->t_grid, m->t_broad, m->t_narrow, m->t_total,
            m->t_M, m->t_N,
            m->kinetic_energy, m->momentum_magnitude);
}

void log_close(FILE *f) {
    if (f) fclose(f);
}
