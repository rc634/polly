#include <iostream>
#include "field.hpp"
#include "grid.hpp"
#include "multigrid.hpp"
#include "problem.hpp"
#include <cmath>


int main() {

    Multigrid mg;

    double old_delta;
    double new_delta;
    double old_integral;

    mg.init();
    mg.hello();

    // v cycle controlls 
    int i_top = mg.i_fine;
    int i_bot = mg.i_coarse;
    int n_trigger = 200;

    // set zeros for the solver to start from
    mg.initial_data();
    mg.save_data("g");

    // give it a tolerance
    mg.up_cycle(10e-12);

    mg.save_data("f");

    // saves with names psi and W
    mg.fine_grid_ptr->save_state();

    return 0;
}

