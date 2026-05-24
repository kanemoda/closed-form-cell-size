/*
 * Phase 8 / Phase 9 Part 1 -- the D2 thesis test on STATIC point sets with a
 * known, tunable correlation dimension, parameterised by dimension and N:
 *
 *     ./phase8_d2 <dim> <N>      (dim = 2 or 3)
 *
 * Does estimating D2 at runtime (Mode B) beat the naive D2:=d (Mode A), and
 * does the benefit SHARPEN with N (the cost-curvature prediction: at small N
 * the ell-independent insertion cost c0 dilutes the penalty; as N grows the
 * ell-dependent cost dominates and the full penalty (~ (1/2) d D2 delta^2)
 * shows)?
 *
 * Configurations (ground-truth D2 exact by self-similar construction):
 *   2D: uniform (D2=2), 1D sine curve (D2~1), 2D Cantor dust keeping k of b^2
 *       sub-cells (D2 = log k/log b), spanning D2 = 0.5 .. 2.0.
 *   3D: uniform (D2=3), a 2D wavy sheet (D2~2), 3D Cantor dust keeping k of
 *       b^3 sub-cells, spanning D2 = 0.5 .. 3.0.
 * Cantor depth ~ ln N / ln k (~1 particle per surviving cell), so the fractal
 * scaling is clean from finest << 2r up through the operative range.
 *
 * Measurement is leak-free (all t_detect via oracle_eval[_split], never the
 * sim_step path). a,b are calibrated once per (dim,N) from a dynamic stable run
 * (a single static snapshot is degenerate: M,S collinear along the ell sweep).
 */

#define _POSIX_C_SOURCE 199309L
#include "rng.h"
#include "config.h"
#include "sim.h"
#include "scenarios.h"
#include "grid.h"
#include "oracle.h"
#include "adapter.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define RADIUS 0.5
#define GAMMA  1.3
#define DOM2D  729.0          /* 2D box (3^6: clean Cantor cells)            */
#define DOM3D  200.0          /* 3D box (matches Phase 6; M(1)=8e6 avoided)  */

typedef struct {
    char    name[24];
    double  analytic_d2;
    double  finest;           /* finest fractal cell (0 for non-fractal)     */
    int     dim, n;
    double  W, H, D;          /* D = depth extent (3D)                       */
    double *px, *py, *pz;     /* pz NULL in 2D                               */
} d2cfg_t;

static int rng_int(rng_state_t *r, int k) { return (int)(rng_next_u64(r) % (uint64_t)k); }

/* ---------------------------------------------------------- generators */

static void gen_uniform(d2cfg_t *c, uint64_t seed) {
    rng_state_t r; rng_seed(&r, seed);
    for (int i = 0; i < c->n; i++) {
        c->px[i] = rng_uniform(&r, 2.0, c->W - 2.0);
        c->py[i] = rng_uniform(&r, 2.0, c->H - 2.0);
        if (c->dim == 3) c->pz[i] = rng_uniform(&r, 2.0, c->D - 2.0);
    }
    c->analytic_d2 = (double)c->dim; c->finest = 0.0;
}

/* 2D: 1D sine curve (D2~1). */
static void gen_curve(d2cfg_t *c, uint64_t seed) {
    rng_state_t r; rng_seed(&r, seed);
    for (int i = 0; i < c->n; i++) {
        double x = rng_uniform(&r, 2.0, c->W - 2.0);
        double y = c->H * 0.5 + 0.30 * c->H * sin(2.0 * 2.0 * M_PI * x / c->W);
        y += rng_uniform(&r, -0.4, 0.4);
        if (y < 2.0) y = 2.0;
        if (y > c->H - 2.0) y = c->H - 2.0;
        c->px[i] = x; c->py[i] = y;
    }
    c->analytic_d2 = 1.0; c->finest = 0.0;
}

/* 3D: a 2D wavy sheet embedded in 3D (D2~2). */
static void gen_sheet3d(d2cfg_t *c, uint64_t seed) {
    rng_state_t r; rng_seed(&r, seed);
    for (int i = 0; i < c->n; i++) {
        double x = rng_uniform(&r, 2.0, c->W - 2.0);
        double y = rng_uniform(&r, 2.0, c->H - 2.0);
        double z = c->D * 0.5 + 0.15 * c->D * sin(2.0 * M_PI * x / c->W) * sin(2.0 * M_PI * y / c->H);
        z += rng_uniform(&r, -0.4, 0.4);
        if (z < 2.0) z = 2.0;
        if (z > c->D - 2.0) z = c->D - 2.0;
        c->px[i] = x; c->py[i] = y; c->pz[i] = z;
    }
    c->analytic_d2 = 2.0; c->finest = 0.0;
}

/* Cantor dust (2D: b^2 sub-cells; 3D: b^3). Keep k spread sub-cells/level;
 * each particle random-walks the kept set down `depth` levels. */
static void gen_cantor(d2cfg_t *c, int b, int k, int depth, uint64_t seed) {
    int nsub = (c->dim == 3) ? b * b * b : b * b;
    int keep[729];
    for (int j = 0; j < k; j++) keep[j] = (j * nsub) / k;
    rng_state_t r; rng_seed(&r, seed);
    for (int i = 0; i < c->n; i++) {
        double x0 = 0, y0 = 0, z0 = 0, sw = c->W, sh = c->H, sd = c->D;
        for (int l = 0; l < depth; l++) {
            int s = keep[rng_int(&r, k)];
            int cx = s % b, cy = (s / b) % b, cz = s / (b * b);
            sw /= b; sh /= b;
            x0 += cx * sw; y0 += cy * sh;
            if (c->dim == 3) { sd /= b; z0 += cz * sd; }
        }
        c->px[i] = x0 + rng_uniform(&r, 0.0, sw);
        c->py[i] = y0 + rng_uniform(&r, 0.0, sh);
        if (c->dim == 3) c->pz[i] = z0 + rng_uniform(&r, 0.0, sd);
    }
    c->analytic_d2 = log((double)k) / log((double)b);
    c->finest = c->W / pow((double)b, (double)depth);
}

/* ------------------------------------------------------- dim-aware probes */

static int num_cells(const d2cfg_t *c, double ell) {
    int gw = (int)ceil(c->W / ell), gh = (int)ceil(c->H / ell);
    if (gw < 1) gw = 1;
    if (gh < 1) gh = 1;
    if (c->dim == 3) {
        int gd = (int)ceil(c->D / ell); if (gd < 1) gd = 1;
        return gw * gh * gd;
    }
    return gw * gh;
}

static long count_S(const d2cfg_t *c, double ell) {
    if (c->dim == 3)
        return grid_count_S_at3d(c->W, c->H, c->D, c->n, c->px, c->py, c->pz, ell);
    return grid_count_S_at(c->W, c->H, c->n, c->px, c->py, ell);
}

/* correlation dimension C(r)~r^D2 (subsampled to <=8000 pts; O(n_s^2)). */
static double correlation_dim(const d2cfg_t *c, double rmin, double rmax) {
    int stride = c->n / 8000; if (stride < 1) stride = 1;
    enum { NB = 24 };
    double lrmin = log(rmin), dl = (log(rmax) - lrmin) / NB;
    long hist[NB]; for (int k = 0; k < NB; k++) hist[k] = 0;
    for (int i = 0; i < c->n; i += stride)
        for (int j = i + stride; j < c->n; j += stride) {
            double dx = c->px[j]-c->px[i], dy = c->py[j]-c->py[i];
            double r2 = dx*dx + dy*dy;
            if (c->dim == 3) { double dz = c->pz[j]-c->pz[i]; r2 += dz*dz; }
            if (r2 <= 0) continue;
            double r = sqrt(r2);
            if (r < rmin || r >= rmax) continue;
            int kk = (int)((log(r) - lrmin) / dl);
            if (kk >= 0 && kk < NB) hist[kk]++;
        }
    long cum = 0; double sx=0, sy=0, sxx=0, sxy=0; int m=0;
    for (int k = 0; k < NB; k++) {
        cum += hist[k];
        if (cum <= 0) continue;
        double rr = exp(lrmin + (k+1)*dl), x = log(rr), y = log((double)cum);
        sx+=x; sy+=y; sxx+=x*x; sxy+=x*y; m++;
    }
    if (m < 2) return 0.0;
    return (m*sxy - sx*sy) / (m*sxx - sx*sx);
}

static double estimator_d2(const d2cfg_t *c, double ell) {
    long Slo = count_S(c, ell / GAMMA), Sm = count_S(c, ell), Shi = count_S(c, ell * GAMMA);
    return adapter_estimate_D2(ell/GAMMA, ell, ell*GAMMA, (double)Slo, (double)Sm, (double)Shi);
}

static double t_at(const d2cfg_t *c, double ell, int K) {
    oracle_point_t p;
    if (c->dim == 3)
        oracle_eval3d(c->n, c->px, c->py, c->pz, RADIUS, c->W, c->H, c->D, ell, K, &p);
    else
        oracle_eval(c->n, c->px, c->py, RADIUS, c->W, c->H, ell, K, &p);
    return p.t_detect * 1e6;
}

static void oracle_best(const d2cfg_t *c, const double *ells, int n_ells, int K,
                        double *best_ell, double *best_t_us) {
    double bsel = 1e300, beval = 0, bell = ells[0];
    for (int e = 0; e < n_ells; e++) {
        oracle_point_t sel, ev;
        if (c->dim == 3)
            oracle_eval_split3d(c->n, c->px, c->py, c->pz, RADIUS, c->W, c->H, c->D, ells[e], K, &sel, &ev);
        else
            oracle_eval_split(c->n, c->px, c->py, RADIUS, c->W, c->H, ells[e], K, &sel, &ev);
        if (sel.t_detect < bsel) { bsel = sel.t_detect; beval = ev.t_detect; bell = ells[e]; }
    }
    *best_ell = bell; *best_t_us = beval * 1e6;
}

/* iterate the closed form to its fixed point; use_est=0 Mode A, 1 Mode B. */
static double mode_fixed_point(const d2cfg_t *c, double a, double b, int use_est,
                               double ell_min, double ell_max) {
    int d = c->dim;
    double ell = pow(c->W * c->H * (d == 3 ? c->D : 1.0) / c->n, 1.0 / d);  /* density baseline */
    ell = adapter_clamp(ell, ell_min, ell_max);
    for (int it = 0; it < 200; it++) {
        int  M = num_cells(c, ell);
        long S = count_S(c, ell);
        double D2 = use_est ? estimator_d2(c, ell) : (double)d;
        D2 = adapter_clamp(D2, 0.2, (double)d);
        double en = adapter_optimal_ell(ell, d, D2, a * M, b * (double)S);
        en = adapter_clamp(en, ell_min, ell_max);
        double step = fabs(log(en / ell));
        ell = en;
        if (step < 1e-5) break;
    }
    return ell;
}

/* --------------------------------------------------------- calibration */

static cost_fit_t calibrate_ab(int dim, int N) {
    config_t cfg; config_default(&cfg);
    cfg.dim = dim; cfg.N = N; cfg.radius = RADIUS; cfg.dt = 1.0/60.0;
    cfg.init_speed = 10.0; cfg.restitution = 0.95; cfg.num_frames = 60;
    cfg.seed = 31; cfg.cell_size = 2.0;
    snprintf(cfg.scenario, sizeof(cfg.scenario), "stable");
    cfg.log_enabled = 0;
    double ladmax;
    if (dim == 3) { cfg.domain_w = cfg.domain_h = cfg.domain_d = DOM3D; ladmax = 50.0; }
    else          { cfg.domain_w = 800; cfg.domain_h = 600; ladmax = 20.0; }
    double ladmin = (dim == 3) ? 1.5 : 1.0;

    double ells[32];
    int n_ells = oracle_make_geo_sweep(ladmin, ladmax, 1.3, ells, 32);
    sim_t s; sim_init(&s, &cfg); scenario_init(&s, SCEN_STABLE);
    int F = 6, K = 5, cap = F * n_ells, ns = 0;
    double *M = malloc(sizeof(double)*cap), *S = malloc(sizeof(double)*cap), *T = malloc(sizeof(double)*cap);
    oracle_point_t *pts = malloc(sizeof(*pts) * n_ells);
    int targets[6]; for (int j = 0; j < F; j++) targets[j] = (j+1)*cfg.num_frames/(F+1);
    int tk = 0;
    for (int t = 0; t < cfg.num_frames && tk < F; t++) {
        sim_step(&s, NULL);
        if (t != targets[tk]) continue;
        if (dim == 3)
            oracle_eval_sweep3d(s.n, s.px, s.py, s.pz, cfg.radius, cfg.domain_w, cfg.domain_h, cfg.domain_d, ells, n_ells, K, pts);
        else
            oracle_eval_sweep(s.n, s.px, s.py, cfg.radius, cfg.domain_w, cfg.domain_h, ells, n_ells, K, pts);
        for (int i = 0; i < n_ells && ns < cap; i++) { M[ns]=pts[i].num_cells; S[ns]=(double)pts[i].S; T[ns]=pts[i].t_detect; ns++; }
        tk++;
    }
    cost_fit_t f = adapter_fit_cost_model(M, S, T, ns);
    free(M); free(S); free(T); free(pts); sim_free(&s);
    return f;
}

/* ----------------------------------------------------------------- main */

static d2cfg_t make(const char *name, int dim, int n) {
    d2cfg_t c; memset(&c, 0, sizeof(c));
    snprintf(c.name, sizeof(c.name), "%s", name);
    c.dim = dim; c.n = n;
    c.W = c.H = (dim == 3) ? DOM3D : DOM2D;
    c.D = (dim == 3) ? DOM3D : 0.0;
    c.px = malloc(sizeof(double)*n); c.py = malloc(sizeof(double)*n);
    c.pz = (dim == 3) ? malloc(sizeof(double)*n) : NULL;
    return c;
}
static int idepth(int N, int k) { int d = (int)(log((double)N)/log((double)k) + 0.5); return d < 2 ? 2 : d; }

int main(int argc, char **argv) {
    int dim = 2, N = 8000, off = 0;
    if (argc > 1) dim = atoi(argv[1]);
    if (argc > 2) N = atoi(argv[2]);
    if (argc > 3) off = atoi(argv[3]);   /* seed offset -> independent realisations */

    d2cfg_t cfgs[8]; int nc = 0;
    if (dim == 2) {
        cfgs[nc]=make("uniform",2,N);     gen_uniform(&cfgs[nc],101+off);                           nc++;
        cfgs[nc]=make("cantor b3k6",2,N); gen_cantor (&cfgs[nc],3,6,idepth(N,6),102+off);           nc++;
        cfgs[nc]=make("cantor b3k4",2,N); gen_cantor (&cfgs[nc],3,4,idepth(N,4),103+off);           nc++;
        cfgs[nc]=make("cantor b3k3",2,N); gen_cantor (&cfgs[nc],3,3,idepth(N,3),104+off);           nc++;
        cfgs[nc]=make("curve",2,N);       gen_curve  (&cfgs[nc],105+off);                           nc++;
        cfgs[nc]=make("cantor b3k2",2,N); gen_cantor (&cfgs[nc],3,2,idepth(N,2),106+off);           nc++;
        cfgs[nc]=make("cantor b9k3",2,N); gen_cantor (&cfgs[nc],9,3,idepth(N,3),107+off);           nc++;
    } else {
        cfgs[nc]=make("uniform",3,N);     gen_uniform(&cfgs[nc],201+off);                           nc++;
        cfgs[nc]=make("cantor b3k9",3,N); gen_cantor (&cfgs[nc],3,9,idepth(N,9),202+off);           nc++;  /* D2=2.0 */
        cfgs[nc]=make("sheet",3,N);       gen_sheet3d(&cfgs[nc],203+off);                           nc++;  /* D2~2  */
        cfgs[nc]=make("cantor b3k5",3,N); gen_cantor (&cfgs[nc],3,5,idepth(N,5),204+off);           nc++;  /* D2=1.46 */
        cfgs[nc]=make("cantor b3k3",3,N); gen_cantor (&cfgs[nc],3,3,idepth(N,3),205+off);           nc++;  /* D2=1.0 */
        cfgs[nc]=make("cantor b3k2",3,N); gen_cantor (&cfgs[nc],3,2,idepth(N,2),206+off);           nc++;  /* D2=0.63 */
        cfgs[nc]=make("cantor b4k2",3,N); gen_cantor (&cfgs[nc],4,2,idepth(N,2),207+off);           nc++;  /* D2=0.5 */
    }

    double ell_min = (dim == 3) ? 1.5 : 1.0;
    /* ladder/clamp max capped above the optima (~8-16) but well below the regime
     * where large cells make S~1e8 (huge pair buffers) at N=100K. */
    double ell_max = (dim == 3) ? 20.0 : 24.0;
    double ells[64];
    /* fine ladder (ratio 1.2) so the discrete-argmin oracle is a tight upper
     * bound vs the modes' continuous ell. */
    int n_ells = oracle_make_geo_sweep(ell_min, ell_max, 1.2, ells, 64);

    fprintf(stderr, "=== Phase 9 Part 1: D2 test, dim=%d N=%d box=%.0f^%d ladder[%.1f,%.1f] ===\n",
            dim, N, (dim==3)?DOM3D:DOM2D, dim, ells[0], ells[n_ells-1]);
    cost_fit_t fit = calibrate_ab(dim, N);
    fprintf(stderr, "calibrated a=%.3e b=%.3e c0=%.3e R^2=%.4f (a/b=%.3f)\n",
            fit.a, fit.b, fit.c0, fit.r2, fit.a/fit.b);
    double a = fit.a, b = fit.b;

    fprintf(stderr, "\n  %-12s | D2_true | corrD2 | estD2@opt(err) | oracle_us oell | ModeA gap(ell) | ModeB gap(ell) | A-B\n", "config");
    fprintf(stderr, "  -------------+---------+--------+----------------+----------------+----------------+----------------+-----\n");
    for (int i = 0; i < nc; i++) {
        double oell, ot; oracle_best(&cfgs[i], ells, n_ells, 3, &oell, &ot);
        double cr_lo = 2.0*RADIUS; if (1.5*cfgs[i].finest > cr_lo) cr_lo = 1.5*cfgs[i].finest;
        double cd  = correlation_dim(&cfgs[i], cr_lo, ((dim==3)?DOM3D:DOM2D)*0.05);
        double ld  = estimator_d2(&cfgs[i], oell);
        double eA  = mode_fixed_point(&cfgs[i], a, b, 0, ell_min, ell_max);
        double eB  = mode_fixed_point(&cfgs[i], a, b, 1, ell_min, ell_max);
        double tA  = t_at(&cfgs[i], eA, 9);
        double tB  = (fabs(eA-eB) < 1e-6) ? tA : t_at(&cfgs[i], eB, 9);
        double gA = tA/ot, gB = tB/ot;
        fprintf(stderr,
            "  %-12s |  %5.3f  | %5.2f  | %5.2f (%+5.2f) | %7.1f %5.2f | %5.2fx(%5.2f) | %5.2fx(%5.2f) | %+.2f\n",
            cfgs[i].name, cfgs[i].analytic_d2, cd, ld, ld-cfgs[i].analytic_d2,
            ot, oell, gA, eA, gB, eB, gA-gB);
        free(cfgs[i].px); free(cfgs[i].py); free(cfgs[i].pz);
    }
    return 0;
}
