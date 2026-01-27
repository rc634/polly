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
    // copy over finest grid 
    fine_grid_ptr = &grids[p.num_grids-1];
    coarse_grid_ptr = &grids[0];

    // set initial data on field
    initial_data();
}

void Multigrid::vcycle() {
    std::cout << "Performing a dummy V-cycle\n";
}

void Multigrid::initial_data() {
    double x0 = 0.2;
    double y0 = 0.3;
    double sig = 0.1;
    fine_grid_ptr->ID_gaussian(x0,y0,sig);

    flush();
}

void Multigrid::flush() {
    // Restriction 

    // loop down over grids to set initial data on all levels
    for (int i = p.num_grids-1; i>= 1; i--)
    {
        restrict_down(grids[i],grids[i-1]);
    }

    // Prolongation 

    // loop up over grids to set initial data on all levels
    for (int i = 0; i < p.num_grids-1; i++)
    {
        prolongate_up(grids[i],grids[i+1]);
    }
}

void Multigrid::v_cycle() {
    // Restriction 

    // loop down over grids to set initial data on all levels
    for (int i = p.num_grids-1; i>= 1; i--)
    {
        grids[i].relax();
        restrict_down(grids[i],grids[i-1]);
    }

    // Prolongation 

    // loop up over grids to set initial data on all levels
    for (int i = 0; i < p.num_grids-1; i++)
    {
        grids[i].relax();
        prolongate_up(grids[i],grids[i+1]);
    }
}

void Multigrid::save_data() {
    fine_grid_ptr->save_data();
    //coarse_grid_ptr->save_data();
}

void Multigrid::hello() {
    std::cout << "Hello from Multigrid with " << p.num_grids << " grids\n";
    // initialize each grid
    for (int i = 0; i < p.num_grids; i++) {
        grids[i].hello();
    }
}

void Multigrid::restrict_down(const Grid &fine, Grid &coarse) {
    for (int j = coarse.f1.jmin; j < coarse.f1.jmax; j++) {
        for (int i = coarse.f1.imin; i < coarse.f1.imax; i++) {
            // fine indices
            int i1 = i*2 - p.ng;
            int i2 = i*2 - p.ng + 1;
            int j1 = j*2 - p.ng;
            int j2 = j*2 - p.ng + 1;

            double val = fine.f1.get_data(i1,j1) + 
                            fine.f1.get_data(i2,j1) + 
                            fine.f1.get_data(i1,j2) + 
                            fine.f1.get_data(i2,j2);

            val *= 0.25;

            coarse.f1.set_data(val,i,j);
        }
    }
}


void Multigrid::prolongate_up(const Grid &coarse, Grid &fine) {
    for (int j = fine.f1.jmin; j < fine.f1.jmax; j++) {
        for (int i = fine.f1.imin; i < fine.f1.imax; i++) {

            // ask for bilinear interpolation 
            double x = fine.f1.get_x(i,j);
            double y = fine.f1.get_y(i,j);
            double val = coarse.f1.bilinear_interp_xy(x,y);
            fine.f1.set_data(val,i,j);
        }
    }
}
