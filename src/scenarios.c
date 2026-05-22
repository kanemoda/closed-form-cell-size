#include "scenarios.h"
#include "rng.h"
#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const char *names_arr[] = {
    "stable", "collapse", "explosion", "vortex",
    "funnel", "pulse", "multicluster", "stream"
};

const char *scenario_name(scenario_id_t id) {
    if ((int)id < 0 || id >= SCEN_COUNT) return "unknown";
    return names_arr[id];
}

scenario_id_t scenario_from_name(const char *name) {
    if (!name || name[0] == '\0') return SCEN_STABLE;
    for (int i = 0; i < SCEN_COUNT; i++) {
        if (strcmp(name, names_arr[i]) == 0) return (scenario_id_t)i;
    }
    fprintf(stderr, "scenario_from_name: unknown '%s' -- defaulting to stable\n", name);
    return SCEN_STABLE;
}

/* ====================================================== placement helpers */

static void place_nonoverlap_rect(sim_t *s, int from, int to,
                                   double xlo, double xhi,
                                   double ylo, double yhi) {
    const double r      = s->cfg.radius;
    const double rr_min = (2.0 * r) * (2.0 * r);
    const int    max_t  = 200000;
    for (int i = from; i < to; i++) {
        int placed = 0;
        for (int t = 0; t < max_t; t++) {
            double x = rng_uniform(&s->rng, xlo, xhi);
            double y = rng_uniform(&s->rng, ylo, yhi);
            int ok = 1;
            for (int j = 0; j < i; j++) {
                double dx = x - s->px[j];
                double dy = y - s->py[j];
                if (dx * dx + dy * dy < rr_min) { ok = 0; break; }
            }
            if (ok) { s->px[i] = x; s->py[i] = y; placed = 1; break; }
        }
        if (!placed) {
            fprintf(stderr,
                "place_nonoverlap_rect: failed at i=%d in [%g,%g]x[%g,%g] (density?)\n",
                i, xlo, xhi, ylo, yhi);
            exit(1);
        }
    }
}

static void place_nonoverlap_disk(sim_t *s, int from, int to,
                                   double cx, double cy, double R) {
    const double r      = s->cfg.radius;
    const double W      = s->cfg.domain_w;
    const double H      = s->cfg.domain_h;
    const double rr_min = (2.0 * r) * (2.0 * r);
    const int    max_t  = 200000;
    for (int i = from; i < to; i++) {
        int placed = 0;
        for (int t = 0; t < max_t; t++) {
            double th  = rng_uniform(&s->rng, 0.0, 2.0 * M_PI);
            double rho = sqrt(rng_uniform01(&s->rng)) * R;
            double x   = cx + rho * cos(th);
            double y   = cy + rho * sin(th);
            if (x < r || x > W - r || y < r || y > H - r) continue;
            int ok = 1;
            for (int j = 0; j < i; j++) {
                double dx = x - s->px[j];
                double dy = y - s->py[j];
                if (dx * dx + dy * dy < rr_min) { ok = 0; break; }
            }
            if (ok) { s->px[i] = x; s->py[i] = y; placed = 1; break; }
        }
        if (!placed) {
            fprintf(stderr,
                "place_nonoverlap_disk: failed at i=%d center=(%g,%g) R=%g\n",
                i, cx, cy, R);
            exit(1);
        }
    }
}

/* ============================================================ scenarios */

/* --- stable: uniform random placement, random velocity directions.
 *     RNG-call order (place_xy then theta, interleaved per particle) is
 *     kept identical to the pre-Phase-3 sim_init_random so legacy tests
 *     remain bit-identical for the same (seed, config). --- */
static void init_stable(sim_t *s) {
    const double r      = s->cfg.radius;
    const double W      = s->cfg.domain_w, H = s->cfg.domain_h;
    const double rr_min = (2.0 * r) * (2.0 * r);
    const int    max_t  = 50000;
    const int    n      = s->n;
    for (int i = 0; i < n; i++) {
        int placed = 0;
        for (int t = 0; t < max_t; t++) {
            double x = rng_uniform(&s->rng, r, W - r);
            double y = rng_uniform(&s->rng, r, H - r);
            int ok = 1;
            for (int j = 0; j < i; j++) {
                double dx = x - s->px[j];
                double dy = y - s->py[j];
                if (dx * dx + dy * dy < rr_min) { ok = 0; break; }
            }
            if (ok) { s->px[i] = x; s->py[i] = y; placed = 1; break; }
        }
        if (!placed) {
            fprintf(stderr, "init_stable: failed to place particle %d (density?)\n", i);
            exit(1);
        }
        double th = rng_uniform(&s->rng, 0.0, 2.0 * M_PI);
        s->vx[i] = s->cfg.init_speed * cos(th);
        s->vy[i] = s->cfg.init_speed * sin(th);
    }
}

/* --- collapse: random placement, inward radial velocity (D2 stays ~d but
 *               the cloud contracts and density rises) --- */
static void init_collapse(sim_t *s) {
    const double r  = s->cfg.radius;
    const double W  = s->cfg.domain_w, H = s->cfg.domain_h;
    const double cx = W * 0.5, cy = H * 0.5;
    place_nonoverlap_rect(s, 0, s->n, r, W - r, r, H - r);
    for (int i = 0; i < s->n; i++) {
        double dx = s->px[i] - cx, dy = s->py[i] - cy;
        double d  = sqrt(dx * dx + dy * dy);
        if (d < 1e-9) {
            double th = rng_uniform(&s->rng, 0.0, 2.0 * M_PI);
            s->vx[i]  = s->cfg.init_speed * cos(th);
            s->vy[i]  = s->cfg.init_speed * sin(th);
        } else {
            s->vx[i] = -s->cfg.init_speed * dx / d;
            s->vy[i] = -s->cfg.init_speed * dy / d;
        }
    }
}

/* --- explosion: cluster in a central disk, outward radial velocity --- */
static void init_explosion(sim_t *s) {
    const double r  = s->cfg.radius;
    const double W  = s->cfg.domain_w, H = s->cfg.domain_h;
    const double cx = W * 0.5, cy = H * 0.5;
    /* Cluster radius: large enough to pack N particles at moderate density.
     * Heuristic R >= r*sqrt(2N) keeps packing fraction below ~0.6.       */
    double R_for_N = r * sqrt(2.0 * (double)s->n) * 1.3;
    double R_max   = 0.45 * fmin(W, H) - r;
    double R       = fmax(fmin(W, H) * 0.1, R_for_N);
    if (R > R_max) R = R_max;
    place_nonoverlap_disk(s, 0, s->n, cx, cy, R);
    for (int i = 0; i < s->n; i++) {
        double dx = s->px[i] - cx, dy = s->py[i] - cy;
        double d  = sqrt(dx * dx + dy * dy);
        if (d < 1e-9) {
            double th = rng_uniform(&s->rng, 0.0, 2.0 * M_PI);
            s->vx[i]  = s->cfg.init_speed * cos(th);
            s->vy[i]  = s->cfg.init_speed * sin(th);
        } else {
            s->vx[i] = s->cfg.init_speed * dx / d;
            s->vy[i] = s->cfg.init_speed * dy / d;
        }
    }
}

/* --- vortex/spiral: 2-arm Archimedean spiral with arm thickness,
 *                    tangential (CCW) initial velocity --- */
static void init_vortex(sim_t *s) {
    const double r  = s->cfg.radius;
    const double W  = s->cfg.domain_w, H = s->cfg.domain_h;
    const double cx = W * 0.5, cy = H * 0.5;
    const int    num_arms      = 2;
    const int    turns         = 3;          /* 3 full revolutions per arm */
    const double max_r         = fmin(W, H) * 0.4 - r;
    const double arm_thickness = max_r * 0.04;
    const double rr_min        = (2.0 * r) * (2.0 * r);
    const int    max_t         = 200000;

    int placed_count = 0;
    for (int i = 0; i < s->n; i++) {
        int arm = i % num_arms;
        int placed = 0;
        for (int t = 0; t < max_t; t++) {
            /* Skip the very center (5%) so the inner arm isn't degenerate. */
            double u       = rng_uniform01(&s->rng);
            double th_par  = (0.05 + 0.95 * u) * 2.0 * M_PI * turns;
            double sp_r    = max_r * (th_par / (2.0 * M_PI * turns));
            double arc_th  = th_par + arm * (2.0 * M_PI / num_arms);
            double jx      = rng_uniform(&s->rng, -arm_thickness, arm_thickness);
            double jy      = rng_uniform(&s->rng, -arm_thickness, arm_thickness);
            double x       = cx + sp_r * cos(arc_th) + jx;
            double y       = cy + sp_r * sin(arc_th) + jy;
            if (x < r || x > W - r || y < r || y > H - r) continue;
            int ok = 1;
            for (int j = 0; j < placed_count; j++) {
                double dx = x - s->px[j], dy = y - s->py[j];
                if (dx * dx + dy * dy < rr_min) { ok = 0; break; }
            }
            if (ok) {
                s->px[placed_count] = x;
                s->py[placed_count] = y;
                /* CCW tangential velocity */
                s->vx[placed_count] = -s->cfg.init_speed * sin(arc_th);
                s->vy[placed_count] =  s->cfg.init_speed * cos(arc_th);
                placed_count++;
                placed = 1;
                break;
            }
        }
        if (!placed) {
            fprintf(stderr, "init_vortex: failed at i=%d (density?)\n", placed_count);
            exit(1);
        }
    }
}

/* --- funnel: sinusoidal attractor along y = y_center + A*sin(B*x).
 *             Spring force pulls every particle toward the curve in y. --- */
static void init_funnel(sim_t *s) {
    const double r = s->cfg.radius;
    const double W = s->cfg.domain_w, H = s->cfg.domain_h;
    place_nonoverlap_rect(s, 0, s->n, r, W - r, r, H - r);
    /* Small random velocities -- the funnel force does most of the work. */
    for (int i = 0; i < s->n; i++) {
        double th = rng_uniform(&s->rng, 0.0, 2.0 * M_PI);
        s->vx[i]  = s->cfg.init_speed * 0.3 * cos(th);
        s->vy[i]  = s->cfg.init_speed * 0.3 * sin(th);
    }
}

static void force_funnel(sim_t *s) {
    const double W  = s->cfg.domain_w, H = s->cfg.domain_h;
    const double dt = s->cfg.dt;
    /* Two wavelengths across the domain (k = 5 keeps oscillation slow
     * relative to the symplectic step: sqrt(5)*dt = 0.037 << 1).         */
    const double k_spring = 5.0;
    const double A        = H * 0.15;
    const double B        = 4.0 * M_PI / W;
    const double y_center = H * 0.5;
    for (int i = 0; i < s->n; i++) {
        double y_target = y_center + A * sin(B * s->px[i]);
        s->vy[i] += -k_spring * (s->py[i] - y_target) * dt;
    }
}

/* --- pulse: a constant attractive central spring.
 *
 *   We interpret the spec's "periodic radial spring" as the natural
 *   periodic oscillation of particles around the centre of attraction
 *   (period pi/sqrt(k) for the half-amplitude / spread). Choosing a
 *   constant k keeps the integrator exactly symplectic and energy
 *   bounded; the *spread* still oscillates with a clean period, which
 *   is what the spec gate (tracking-lag stress) actually needs. --- */
static void init_pulse(sim_t *s) {
    const double r = s->cfg.radius;
    const double W = s->cfg.domain_w, H = s->cfg.domain_h;
    place_nonoverlap_rect(s, 0, s->n, r, W - r, r, H - r);
    for (int i = 0; i < s->n; i++) {
        double th = rng_uniform(&s->rng, 0.0, 2.0 * M_PI);
        s->vx[i]  = s->cfg.init_speed * 0.05 * cos(th);
        s->vy[i]  = s->cfg.init_speed * 0.05 * sin(th);
    }
}

static void force_pulse(sim_t *s) {
    const double W  = s->cfg.domain_w, H = s->cfg.domain_h;
    const double dt = s->cfg.dt;
    const double cx = W * 0.5, cy = H * 0.5;
    /* omega_natural = sqrt(k) = 2 rad/s. spread-period = pi/omega ≈ 1.57 s. */
    const double k_spring = 4.0;
    for (int i = 0; i < s->n; i++) {
        s->vx[i] += -k_spring * (s->px[i] - cx) * dt;
        s->vy[i] += -k_spring * (s->py[i] - cy) * dt;
    }
}

/* --- multicluster: K clusters evenly arranged on a ring, each cluster
 *                   travelling inward toward the global centre --- */
static void init_multicluster(sim_t *s) {
    const double W = s->cfg.domain_w, H = s->cfg.domain_h;
    const int    K = 4;
    const double ring_r    = fmin(W, H) * 0.25;
    double       cluster_R = ring_r * 0.4;
    double       cap       = fmin(W, H) * 0.15;
    if (cluster_R > cap) cluster_R = cap;
    const int per = s->n / K;
    int idx = 0;
    for (int k = 0; k < K; k++) {
        int end = (k == K - 1) ? s->n : (idx + per);
        double th = 2.0 * M_PI * (double)k / (double)K;
        double cx = W * 0.5 + ring_r * cos(th);
        double cy = H * 0.5 + ring_r * sin(th);
        place_nonoverlap_disk(s, idx, end, cx, cy, cluster_R);
        double dx_to_center = W * 0.5 - cx;
        double dy_to_center = H * 0.5 - cy;
        double d            = sqrt(dx_to_center * dx_to_center + dy_to_center * dy_to_center);
        for (int i = idx; i < end; i++) {
            s->vx[i] = s->cfg.init_speed * dx_to_center / d;
            s->vy[i] = s->cfg.init_speed * dy_to_center / d;
        }
        idx = end;
    }
}

/* --- stream: particles in a thin horizontal band, x-velocity = +init_speed,
 *             with a y-confining spring keeping them in the channel --- */
static void init_stream(sim_t *s) {
    const double r = s->cfg.radius;
    const double W = s->cfg.domain_w, H = s->cfg.domain_h;
    const double channel_half = H * 0.1;
    const double y_center     = H * 0.5;
    place_nonoverlap_rect(s, 0, s->n, r, W - r,
                          y_center - channel_half, y_center + channel_half);
    for (int i = 0; i < s->n; i++) {
        s->vx[i] = s->cfg.init_speed;
        s->vy[i] = 0.0;
    }
}

static void force_stream(sim_t *s) {
    const double H        = s->cfg.domain_h;
    const double dt       = s->cfg.dt;
    const double y_center = H * 0.5;
    const double k_y      = 3.0;     /* gentle confining spring; sqrt(3)*dt ≈ 0.029 */
    for (int i = 0; i < s->n; i++) {
        s->vy[i] += -k_y * (s->py[i] - y_center) * dt;
    }
}

/* ================================================================ vtable */
typedef struct {
    void (*init)(sim_t *);
    void (*apply_forces)(sim_t *);   /* may be NULL */
} scenario_vt;

static const scenario_vt vtables[SCEN_COUNT] = {
    [SCEN_STABLE]       = { init_stable,       NULL          },
    [SCEN_COLLAPSE]     = { init_collapse,     NULL          },
    [SCEN_EXPLOSION]    = { init_explosion,    NULL          },
    [SCEN_VORTEX]       = { init_vortex,       NULL          },
    [SCEN_FUNNEL]       = { init_funnel,       force_funnel  },
    [SCEN_PULSE]        = { init_pulse,        force_pulse   },
    [SCEN_MULTICLUSTER] = { init_multicluster, NULL          },
    [SCEN_STREAM]       = { init_stream,       force_stream  },
};

int scenario_has_forces(scenario_id_t id) {
    if ((int)id < 0 || id >= SCEN_COUNT) return 0;
    return vtables[id].apply_forces != NULL;
}

void scenario_init(sim_t *s, scenario_id_t id) {
    if ((int)id < 0 || id >= SCEN_COUNT) {
        fprintf(stderr, "scenario_init: bad id %d\n", (int)id);
        exit(1);
    }
    s->scenario_id = (int)id;
    s->frame_index = 0;
    vtables[id].init(s);
}

void scenario_apply_forces(sim_t *s) {
    int id = s->scenario_id;
    if (id < 0 || id >= SCEN_COUNT) return;
    if (vtables[id].apply_forces) vtables[id].apply_forces(s);
}
