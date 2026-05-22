#ifndef LOG_H
#define LOG_H

#include <stdio.h>
#include "sim.h"

/*
 * Per-frame CSV log. One row per frame, columns:
 *   frame, cell_size, N, num_cells, S, candidate_pairs, collisions,
 *   t_grid, t_broad, t_narrow, t_total, kinetic_energy, momentum_magnitude
 *
 * Floating-point columns are written with %.17g (round-trip-safe IEEE-754),
 * timer columns with %.9g (nanosecond precision is enough for analysis).
 * The deterministic columns (everything except the four timers) are
 * bit-identical across runs of the same (config, seed).
 */
FILE *log_open(const char *path);
void  log_write_header(FILE *f);
void  log_write_row   (FILE *f, const sim_metrics_t *m);
void  log_close       (FILE *f);

#endif
