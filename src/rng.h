#ifndef RNG_H
#define RNG_H

#include <stdint.h>

/*
 * Named, seeded RNG: xoshiro256++ (Blackman & Vigna, 2018).
 * State seeded via SplitMix64 from a single uint64_t.
 * Same seed -> bit-identical stream, on any machine that does
 * IEEE-754 double arithmetic.
 */

typedef struct {
    uint64_t s[4];
} rng_state_t;

void     rng_seed(rng_state_t *st, uint64_t seed);
uint64_t rng_next_u64(rng_state_t *st);
double   rng_uniform01(rng_state_t *st);                          /* [0, 1) */
double   rng_uniform(rng_state_t *st, double lo, double hi);      /* [lo, hi) */

#endif
