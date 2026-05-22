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
} config_t;

void config_default(config_t *cfg);
int  config_load   (const char *path, config_t *cfg);
void config_print  (const config_t *cfg);

#endif
