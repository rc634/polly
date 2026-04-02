#include <iostream>
#include "field.hpp"
#include "grid.hpp"
#include "multigrid.hpp"
#include "problem.hpp"
#include "params.hpp"
#include <cmath>


int main() {
    
    Params p; 
    Multigrid mg;

    double old_delta;
    double new_delta;
    double old_integral;

    // automatically handles wether analytic source or not
    mg.init(); 

    // initial check
    mg.hello();

    // set zeros for the solver to start from
    mg.initial_data();
    mg.save_data("g");

    // give it a tolerance
    mg.up_cycle(p.tol);

    mg.save_data("f");

    // saves with names psi and W
    mg.fine_grid_ptr->save_state();

    return 0;
}

