#ifndef CONFIG_H
#define CONFIG_H

#include <stdint.h>

typedef struct {
    int      N;             /* particle count                       */
    double   domain_w;      /* domain width (x extent)              */
    double   domain_h;      /* domain height (y extent)             */
    double   radius;        /* particle radius (equal for all)      */
    double   restitution;   /* pair collision restitution (0..1)    */
    double   dt;            /* fixed integrator timestep            */
    double   cell_size;     /* uniform-grid cell edge length (>=2r) */
    int      num_frames;    /* total frames to simulate             */
    uint64_t seed;          /* RNG seed                             */
    double   init_speed;    /* initial velocity magnitude           */

    /* Phase 2: per-frame CSV logging */
    int      log_enabled;            /* 0 = disabled, 1 = enabled    */
    char     log_csv_path[256];      /* output CSV path; "" disables */

    /* Phase 3: scenario selection + optional snapshot dump.
     * scenario: one of "stable", "collapse", "explosion", "vortex",
     *           "funnel", "pulse", "multicluster", "stream".
     * snapshot_csv_path: if non-empty, main.c dumps particle positions
     *           at 5 evenly spaced frames to this CSV (frame,id,x,y,vx,vy). */
    char     scenario[32];
    char     snapshot_csv_path[256];

    /* Phase 4: per-frame oracle (non-perturbing snapshot sweep).
     * oracle_csv_path: empty disables the per-frame oracle entirely;
     *           otherwise main.c writes (frame, ell, S, num_cells,
     *           candidate_pairs, t_grid, t_broad, t_narrow, t_detect)
     *           for every (frame, ell) pair to this CSV.
     * oracle_ell_min/max: bounds of the geometric sweep. <=0 -> defaults
     *           (2r and min(W,H)/2 respectively).
     * oracle_gamma:  geometric ratio between successive ells. <=1 -> 1.25.
     * oracle_K:      median-over-K timing repeats per ell.  <1   -> 5.    */
    char     oracle_csv_path[256];
    double   oracle_ell_min;
    double   oracle_ell_max;
    double   oracle_gamma;
    int      oracle_K;
} config_t;

void config_default(config_t *cfg);
int  config_load   (const char *path, config_t *cfg);
void config_print  (const config_t *cfg);

#endif
