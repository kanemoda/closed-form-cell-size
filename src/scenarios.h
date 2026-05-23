#ifndef SCENARIOS_H
#define SCENARIOS_H

#include "sim.h"

/*
 * Phase 3 scenarios (spec section 9.1).
 *
 * Each scenario is a deterministic function of (seed, config): it sets the
 * initial positions and velocities, and -- where the spec calls for it --
 * supplies a per-frame external force field that the semi-implicit-Euler
 * force pass in sim_step applies before the position update.
 *
 * Force-field design: F enters via v <- v + (F/m)*dt during sim_step's
 * force pass. With m = 1, this is just v += F*dt. Per-scenario forces are
 * deliberately bounded so the symplectic integrator keeps total energy
 * (KE + PE) bounded over long runs.
 */
typedef enum {
    SCEN_STABLE       = 0,   /* uniform random, no force                      */
    SCEN_COLLAPSE     = 1,   /* random placement, inward radial velocity      */
    SCEN_EXPLOSION    = 2,   /* central cluster, outward radial velocity      */
    SCEN_VORTEX       = 3,   /* 2-arm Archimedean spiral, tangential v        */
    SCEN_FUNNEL       = 4,   /* sinusoidal y-attractor (filament formation)   */
    SCEN_PULSE        = 5,   /* central attractive spring (radial oscillation)*/
    SCEN_MULTICLUSTER = 6,   /* K clusters converging toward common center    */
    SCEN_STREAM       = 7,   /* horizontal channel with y-confining force     */
    SCEN_SETTLING     = 8,   /* gravity-driven settling into a floor layer    */
                             /* (D2 -> d-1 drop test; spec's 3D scenario).    */
    SCEN_COUNT
} scenario_id_t;

/* Spec §9.1 scenario groups. 2D runs the first 8; 3D runs
 * {stable, collapse, explosion, pulse, settling}. scenario_init /
 * scenario_apply_forces dispatch on the sim's cfg.dim, so the same enum value
 * selects the 2D or 3D variant. The 3D-unsupported scenarios (vortex, funnel,
 * multicluster, stream) have no 3D initialiser. */

scenario_id_t scenario_from_name(const char *name);
const char   *scenario_name      (scenario_id_t id);

/* Initialise particle positions/velocities for the given scenario.
 * Uses s->rng (which sim_init seeds from cfg.seed); deterministic.
 * Sets s->scenario_id = id and resets s->frame_index = 0.            */
void scenario_init(sim_t *s, scenario_id_t id);

/* Apply scenario-specific external forces to v. Dispatches on s->scenario_id.
 * No-op for force-free scenarios (STABLE, COLLAPSE, EXPLOSION, VORTEX,
 * MULTICLUSTER). Called from sim_step's symplectic force pass.        */
void scenario_apply_forces(sim_t *s);

/* True if scenario has a non-NULL force field. */
int scenario_has_forces(scenario_id_t id);

#endif
