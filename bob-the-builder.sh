CXX=g++
CXX_OMP=g++-15
FLAGS="-std=c++17 -O2"
FLAGS_OMP="-std=c++17 -O2 -fopenmp"

# ── serial shared objects ──────────────────────────────────────
$CXX $FLAGS -c multigrid.cpp -o multigrid.o
$CXX $FLAGS -c grid.cpp      -o grid.o
$CXX $FLAGS -c field.cpp     -o field.o
$CXX $FLAGS -c source.cpp    -o source.o
SHARED="source.o multigrid.o grid.o field.o"

# polly — single run, serial
$CXX $FLAGS -c main.cpp -o main.o
$CXX $FLAGS main.o $SHARED -o polly
rm main.o

# polly_batch — parameter sweep, serial
$CXX $FLAGS -c main_batch.cpp -o main_batch.o
$CXX $FLAGS main_batch.o $SHARED -o polly_batch
rm main_batch.o

rm $SHARED

# ── OpenMP shared objects ──────────────────────────────────────
$CXX_OMP $FLAGS_OMP -c multigrid.cpp -o multigrid.o
$CXX_OMP $FLAGS_OMP -c grid.cpp      -o grid.o
$CXX_OMP $FLAGS_OMP -c field.cpp     -o field.o
$CXX_OMP $FLAGS_OMP -c source.cpp    -o source.o
SHARED_OMP="source.o multigrid.o grid.o field.o"

# polly_omp — single run, OpenMP
$CXX_OMP $FLAGS_OMP -c main.cpp -o main.o
$CXX_OMP $FLAGS_OMP main.o $SHARED_OMP -o polly_omp
rm main.o

# polly_batch_omp — parameter sweep, OpenMP
$CXX_OMP $FLAGS_OMP -c main_batch.cpp -o main_batch.o
$CXX_OMP $FLAGS_OMP main_batch.o $SHARED_OMP -o polly_batch_omp
rm main_batch.o

rm $SHARED_OMP
