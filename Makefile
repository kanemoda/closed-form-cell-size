CC      ?= gcc
CFLAGS  := -O2 -march=native -Wall -Wextra -Wno-unused-parameter -std=c11 -Isrc
DEPFLAGS := -MMD -MP
LDFLAGS := -lm

CORE_SRC := src/rng.c src/json.c src/config.c src/grid.c src/sim.c src/log.c src/scenarios.c src/oracle.c src/baselines.c src/adapter.c
CORE_OBJ := $(CORE_SRC:.c=.o)

BINS  := simulator test_grid_vs_bruteforce test_physics test_phase2 test_scenarios test_phase4 test_phase5
TOOLS := compare_phase5

.PHONY: all tools test clean
all: $(BINS)
tools: $(TOOLS)

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

test_phase4: $(CORE_OBJ) tests/test_phase4.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test_phase5: $(CORE_OBJ) tests/test_phase5.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

compare_phase5: $(CORE_OBJ) tools/compare_phase5.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# -MMD -MP emits a .d alongside each .o listing the headers it includes, so a
# header edit forces the dependent .o (and binary) to rebuild. Without this a
# stale .o linked against new headers can produce a phantom gate failure.
%.o: %.c
	$(CC) $(CFLAGS) $(DEPFLAGS) -c -o $@ $<

test: $(BINS)
	./test_grid_vs_bruteforce
	./test_physics
	./test_phase2
	./test_scenarios
	./test_phase4
	./test_phase5

clean:
	rm -f src/*.o tests/*.o tools/*.o src/*.d tests/*.d tools/*.d $(BINS) $(TOOLS)

-include $(wildcard src/*.d tests/*.d tools/*.d)
