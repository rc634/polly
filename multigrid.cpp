#include "multigrid.hpp"

Multigrid::Multigrid() {
    // resize vector
    grids.resize(p.num_grids);
    i_coarse = 0;
    i_fine = p.num_grids-1;
}

void Multigrid::init() {
    // initialize each grid
    for (int i = 0; i < p.num_grids; i++) {
        grids[i].init(i);
    }
    // copy over finest grid 
    fine_grid_ptr = &grids[p.num_grids-1];
    coarse_grid_ptr = &grids[0];
}

void Multigrid::fill_all_ghosts() {
    for (int i = 0; i < p.num_grids; i++)
    {
        grids[i].fill_all_ghosts();
    }
}

void Multigrid::initial_data() {
    std::cout << "Zeroing initial data !\n";
    fine_grid_ptr->ID_zeros();
    std::cout << "Integral of solution : " << fine_grid_ptr->field_integral() << "\n";

    // fill_all_ghosts();
}

void Multigrid::analytic_solution() {
    std::cout << "Setting Poisson Disc initial data !\n";
    // loop all grids
    for (size_t i = i_coarse; i <= i_fine; i++)
    {
        grids[i].ID_pdisc(0.1,100); // R_0, rho_0, radius, density
    }
     
    std::cout << "Integral of solution : " << fine_grid_ptr->field_integral() << "\n";
}

// for testing only, can cause large errors
// restricts data from fine to coarse then prolongates data back to fine
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

/////////////////////////////
// V CYCLE 

void Multigrid::v_cycle(int i_top, int i_bot) {
    // Restriction 
    if (i_top > i_fine || i_bot < i_coarse) {
        std::cout << "Warning! incorrect i_top or i_bot given to Multigrid::v_cycle!\n";
        return;
    }

    // loop down over grids to set initial data on all levels
    for (int i = i_top; i > i_bot; i--)
    {
        //std::cout << "Relax grid " << i << "\n"; 
        for (int q = 0; q<p.iter; q++)
        {
            grids[i].relax();
        }
        //std::cout << "Restrict " << i << " to " << i-1 <<"\n";
        restrict_down(grids[i],grids[i-1]);
    }

    // Prolongation 

    // loop up over grids to set initial data on all levels
    for (int i = i_bot; i < i_top; i++)
    {
        //std::cout << "Relax grid " << i << "\n"; 
        for (int q = 0; q<p.iter; q++)
        {
            grids[i].relax();
        }
        //std::cout << "Prolong " << i << " to " << i+1 <<"\n";
        prolongate_up(grids[i],grids[i+1]);
    }

    //std::cout << "Fine level relax\n"; 
    // final smoothing
    for (int q = 0; q<p.iter; q++)
    {
        fine_grid_ptr->relax();
    }
}

/////////////////////////////
// W CYCLE 

void Multigrid::w_cycle() {
    int n = 1;
    for (size_t i = i_coarse; i < i_fine; i++) {
        // q=n repitions 
        for (size_t q = 0; q < n; q++) {
            v_cycle(i_fine, i);
        }
        n*=2;
    }
}

void Multigrid::refine() {
    // Restriction 

    fine_grid_ptr->relax();
}

void Multigrid::save_data(const std::string& filename) {
    fine_grid_ptr->save_data(filename);
    //coarse_grid_ptr->save_data(filename);
}

void Multigrid::hello() {
    std::cout << "Hello from Multigrid with " << p.num_grids << " grids\n";
    // initialize each grid
    for (int i = 0; i < p.num_grids; i++) {
        grids[i].hello();
    }
}

void Multigrid::restrict_down(const Grid &fine, Grid &coarse) {
    for (int j = coarse.jmin; j < coarse.jmax; j++) {
        for (int i = coarse.imin; i < coarse.imax; i++) {
            // fine indices
            int i1 = i*2 - p.ng;
            int i2 = i*2 - p.ng + 1;
            int j1 = j*2 - p.ng;
            int j2 = j*2 - p.ng + 1;

            double val = fine.W.get_data(i1,j1) + 
                            fine.W.get_data(i2,j1) + 
                            fine.W.get_data(i1,j2) + 
                            fine.W.get_data(i2,j2);

            val *= 0.25;

            coarse.W.set_data(val,i,j);
        }
    }
}


void Multigrid::prolongate_up(const Grid &coarse, Grid &fine) {
    for (int j = fine.jmin; j < fine.jmax; j++) {
        for (int i = fine.imin; i < fine.imax; i++) {

            // ask for bilinear interpolation 
            double x = fine.W.get_x(i,j);
            double y = fine.W.get_y(i,j);
            double val = coarse.W.bilinear_interp_xy(x,y);
            fine.W.set_data(val,i,j);
        }
    }
}
