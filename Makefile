CC      ?= gcc
CFLAGS  := -O2 -march=native -Wall -Wextra -Wno-unused-parameter -std=c11 -Isrc
LDFLAGS := -lm

CORE_SRC := src/rng.c src/json.c src/config.c src/grid.c src/sim.c src/log.c src/scenarios.c
CORE_OBJ := $(CORE_SRC:.c=.o)

BINS := simulator test_grid_vs_bruteforce test_physics test_phase2 test_scenarios

.PHONY: all test clean
all: $(BINS)

simulator: $(CORE_OBJ) src/main.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test_grid_vs_bruteforce: $(CORE_OBJ) tests/test_grid_vs_bruteforce.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test_physics: $(CORE_OBJ) tests/test_physics.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test_phase2: $(CORE_OBJ) tests/test_phase2.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test_scenarios: $(CORE_OBJ) tests/test_scenarios.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

test: $(BINS)
	./test_grid_vs_bruteforce
	./test_physics
	./test_phase2
	./test_scenarios

clean:
	rm -f src/*.o tests/*.o $(BINS)
