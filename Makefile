CXX     = g++
CXX_OMP = g++

# -flto lets the accessors in field.cpp / source.cpp inline into the relax
# loop in grid.cpp - without it every get_data / d1x / get_src_* is an
# out-of-line call per cell per step and nothing vectorises
LTO = -flto

FLAGS     = -std=c++17 -O2 $(LTO)
FLAGS_OMP = -std=c++17 -O2 $(LTO) -Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include -L/usr/local/opt/libomp/lib -lomp

SRCS = field.cpp grid.cpp multigrid.cpp source.cpp

# ── serial ────────────────────────────────────────────────────────────────────

polly: main.cpp $(SRCS)
	$(CXX) $(FLAGS) $(SRCS) main.cpp -o $@

polly_batch: main_batch.cpp $(SRCS)
	$(CXX) $(FLAGS) $(SRCS) main_batch.cpp -o $@

all: serial omp

serial: polly polly_batch

# ── OpenMP ────────────────────────────────────────────────────────────────────

polly_omp: main.cpp $(SRCS)
	$(CXX_OMP) $(FLAGS_OMP) $(SRCS) main.cpp -o $@

polly_batch_omp: main_batch.cpp $(SRCS)
	$(CXX_OMP) $(FLAGS_OMP) $(SRCS) main_batch.cpp -o $@

omp: polly_omp polly_batch_omp

# ── clean ─────────────────────────────────────────────────────────────────────

clean:
	rm -f polly polly_batch polly_omp polly_batch_omp

.PHONY: all serial omp clean

