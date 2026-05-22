#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "config.h"
#include "sim.h"

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s <config.json>\n", argv[0]);
        return 1;
    }
    config_t cfg;
    if (config_load(argv[1], &cfg) != 0) return 1;
    config_print(&cfg);

    sim_t s;
    sim_init(&s, &cfg);
    sim_init_random(&s);

    double ke0 = sim_total_kinetic_energy(&s);
    double pxt0, pyt0;
    sim_total_momentum(&s, &pxt0, &pyt0);

    double t0 = now_seconds();
    long total_pairs = 0;
    for (int t = 0; t < cfg.num_frames; t++) {
        sim_step(&s);
        total_pairs += s.last_pair_count;
    }
    double t1 = now_seconds();

    double ke1 = sim_total_kinetic_energy(&s);
    double pxt1, pyt1;
    sim_total_momentum(&s, &pxt1, &pyt1);

    fprintf(stderr,
        "frames=%d  wall=%.3fs  ms/frame=%.4f  total_pairs=%ld  M=%d  S=%ld\n",
        cfg.num_frames, t1 - t0, (t1 - t0) * 1000.0 / cfg.num_frames,
        total_pairs, grid_M(&s.grid), grid_S(&s.grid));
    fprintf(stderr,
        "KE: initial=%.6g final=%.6g  ratio=%.9f\n", ke0, ke1, ke1 / ke0);
    fprintf(stderr,
        "P : initial=(%.6g, %.6g)  final=(%.6g, %.6g)\n", pxt0, pyt0, pxt1, pyt1);

    sim_free(&s);
    return 0;
}
