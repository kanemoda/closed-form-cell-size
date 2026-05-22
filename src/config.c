#include "config.h"
#include "json.h"
#include <stdio.h>
#include <string.h>

void config_default(config_t *cfg) {
    cfg->N           = 1000;
    cfg->domain_w    = 800.0;
    cfg->domain_h    = 600.0;
    cfg->radius      = 0.5;
    cfg->restitution = 0.95;
    cfg->dt          = 1.0 / 60.0;
    cfg->cell_size   = 2.0;     /* = 4r at default r -- safely above the 2r floor */
    cfg->num_frames  = 100;
    cfg->seed        = 1;
    cfg->init_speed  = 50.0;

    cfg->log_enabled = 0;
    cfg->log_csv_path[0] = '\0';

    /* Phase 3 */
    strncpy(cfg->scenario, "stable", sizeof(cfg->scenario) - 1);
    cfg->scenario[sizeof(cfg->scenario) - 1] = '\0';
    cfg->snapshot_csv_path[0] = '\0';

    /* Phase 4 */
    cfg->oracle_csv_path[0] = '\0';
    cfg->oracle_ell_min     = 0.0;
    cfg->oracle_ell_max     = 0.0;
    cfg->oracle_gamma       = 0.0;
    cfg->oracle_K           = 0;
}

int config_load(const char *path, config_t *cfg) {
    config_default(cfg);
    json_object_t obj;
    if (json_parse_file(path, &obj) != 0) {
        fprintf(stderr, "config_load: failed to parse '%s'\n", path);
        return -1;
    }
    json_get_int   (&obj, "N",           &cfg->N);
    json_get_double(&obj, "domain_w",    &cfg->domain_w);
    json_get_double(&obj, "domain_h",    &cfg->domain_h);
    json_get_double(&obj, "radius",      &cfg->radius);
    json_get_double(&obj, "restitution", &cfg->restitution);
    json_get_double(&obj, "dt",          &cfg->dt);
    json_get_double(&obj, "cell_size",   &cfg->cell_size);
    json_get_int   (&obj, "num_frames",  &cfg->num_frames);
    json_get_u64   (&obj, "seed",        &cfg->seed);
    json_get_double(&obj, "init_speed",  &cfg->init_speed);

    json_get_bool  (&obj, "log_enabled", &cfg->log_enabled);
    json_get_string(&obj, "log_csv_path", cfg->log_csv_path,
                    (int)sizeof(cfg->log_csv_path));

    json_get_string(&obj, "scenario",          cfg->scenario,
                    (int)sizeof(cfg->scenario));
    json_get_string(&obj, "snapshot_csv_path", cfg->snapshot_csv_path,
                    (int)sizeof(cfg->snapshot_csv_path));

    json_get_string(&obj, "oracle_csv_path",   cfg->oracle_csv_path,
                    (int)sizeof(cfg->oracle_csv_path));
    json_get_double(&obj, "oracle_ell_min",   &cfg->oracle_ell_min);
    json_get_double(&obj, "oracle_ell_max",   &cfg->oracle_ell_max);
    json_get_double(&obj, "oracle_gamma",     &cfg->oracle_gamma);
    json_get_int   (&obj, "oracle_K",         &cfg->oracle_K);
    return 0;
}

void config_print(const config_t *cfg) {
    fprintf(stderr,
            "config: N=%d domain=%gx%g r=%g e=%g dt=%g cell=%g frames=%d seed=%llu init_speed=%g\n",
            cfg->N, cfg->domain_w, cfg->domain_h, cfg->radius, cfg->restitution,
            cfg->dt, cfg->cell_size, cfg->num_frames,
            (unsigned long long)cfg->seed, cfg->init_speed);
    fprintf(stderr,
            "        log_enabled=%d log_csv_path='%s'\n",
            cfg->log_enabled, cfg->log_csv_path);
    fprintf(stderr,
            "        scenario='%s' snapshot_csv_path='%s'\n",
            cfg->scenario, cfg->snapshot_csv_path);
}
