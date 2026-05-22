#include "rng.h"

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t splitmix64(uint64_t *st) {
    uint64_t z = (*st += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

void rng_seed(rng_state_t *st, uint64_t seed) {
    uint64_t sm = seed;
    for (int i = 0; i < 4; i++) {
        st->s[i] = splitmix64(&sm);
    }
}

uint64_t rng_next_u64(rng_state_t *st) {
    const uint64_t result = rotl(st->s[0] + st->s[3], 23) + st->s[0];
    const uint64_t t = st->s[1] << 17;

    st->s[2] ^= st->s[0];
    st->s[3] ^= st->s[1];
    st->s[1] ^= st->s[2];
    st->s[0] ^= st->s[3];
    st->s[2] ^= t;
    st->s[3] = rotl(st->s[3], 45);

    return result;
}

double rng_uniform01(rng_state_t *st) {
    /* Use the top 53 bits as a uniform double in [0, 1). */
    return (rng_next_u64(st) >> 11) * (1.0 / 9007199254740992.0);
}

double rng_uniform(rng_state_t *st, double lo, double hi) {
    return lo + (hi - lo) * rng_uniform01(st);
}
