#include "multigrid.hpp"

Multigrid::Multigrid() {
    // resize vector
    grids.resize(p.num_grids);
}

void Multigrid::init() {
    // initialize each grid
    for (int i = 0; i < p.num_grids; i++) {
        grids[i].init(i);
    }
}

void Multigrid::vcycle() {
    std::cout << "Performing a dummy V-cycle\n";
}

void Multigrid::save_data() {
    grids[p.num_grids-1].save_data();
}

void Multigrid::hello() {
    std::cout << "Hello from Multigrid with " << p.num_grids << " grids\n";
    // initialize each grid
    for (int i = 0; i < p.num_grids; i++) {
        grids[i].hello();
    }
}
