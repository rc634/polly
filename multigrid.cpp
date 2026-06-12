#include "multigrid.hpp"
#include <cmath>
#include <iomanip>

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

    if (!Source::ANALYTIC) {
        // overwrite all levels with tabulated source from file
        Grid& fg = grids[i_fine];
        fg.m_source.init_from_file("data/source-W.dat", "data/source-psi.dat",
                                    fg.nxg, fg.nyg, fg.ng);
        for (int i = i_fine; i > (int)i_coarse; i--)
            restrict_source_down(grids[i], grids[i-1]);
    }
}

void Multigrid::restrict_source_down(const Grid& fine, Grid& coarse) {
    for (int j = coarse.jmin; j < coarse.jmax; j++) {
        for (int i = coarse.imin; i < coarse.imax; i++) {
            int i1 = i*2 - p.ng, i2 = i*2 - p.ng + 1;
            int j1 = j*2 - p.ng, j2 = j*2 - p.ng + 1;
            coarse.m_source.set_v(i, j,
                0.25 * (fine.m_source.get_v(i1,j1) + fine.m_source.get_v(i2,j1) +
                        fine.m_source.get_v(i1,j2) + fine.m_source.get_v(i2,j2)));
            coarse.m_source.set_rho(i, j,
                0.25 * (fine.m_source.get_rho(i1,j1) + fine.m_source.get_rho(i2,j1) +
                        fine.m_source.get_rho(i1,j2) + fine.m_source.get_rho(i2,j2)));
        }
    }
}

void Multigrid::fill_all_ghosts() {
    for (int i = 0; i < p.num_grids; i++)
    {
        grids[i].fill_all_ghosts();
    }
}


void Multigrid::initial_data() {
    std::cout << "Zeroing initial data !\n";
    for (int i = 0; i < p.num_grids; i++)
    {
        grids[i].ID_zeros();
    }
    // std::cout << "Integral of solution : " << fine_grid_ptr->field_integral() << "\n";

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

// /////////////////////////////
// // UP CYCLE 

// // solves from bottom grid up
// // satisfies accuracy criterion on each grid advancement

// void Multigrid::up_cycle(const double epsilon) {
//     // output precision 
//     std::cout << std::setprecision(8);
//     std::cout << "* -------------------*\n";
//     std::cout << "| STARTING UP CYCLE! |\n";
//     std::cout << "* -------------------*\n";

//     // flush values down so that all levels are initialised
//     flush();
//     double diff = 100.; // init large to be safe
//     int count = 0;
//     double sor = 1.; // over relax param 1-4 in theory
//     double sor_max = 1.5; // in practice max to this!
//     double global_alpha = 1.; // global source blend knob
//     // loop up through all grids
//     for (size_t i = i_coarse; i <= i_fine; i++)
//     {
//         // lazy boi
//         grids[i].fill_all_ghosts();
//         count = 0;

//         std::cout << "Entering Grid Number " << i << std::endl;
//         while (true) {
//             grids[i].m_source.m_alpha = global_alpha;
//             grids[i].m_SOR = sor;
//             for (int q = 0; q<p.iter; q++)
//             {
//                 grids[i].relax();
//             }
//             // sor = sor + 0.2 * (sor_max-sor); // this can bug out, but can be fast
//             diff = abs(grids[i].m_delta);
//             if (diff < epsilon) break;
//             std::cout << " ["<<i<<"/"<<i_fine<<"] - " << count * p.iter << " - diff : " << diff << "\n";
//             std::cout << " - ~ CFL : " << p.CFL*sor << "\n";
//             std::cout << " - ~ int W = " << grids[i].int_W() 
//                       << ", int psi = " << grids[i].int_psi() << "\n";
//             std::cout << " - ~ alpha-blend : " << grids[i].m_source.m_alpha << "\n";
//             std::cout << " * - ~ - ~ - \n";
//             count += 1;
//         }

//         // this creates noise! need to smooth before cranking over relaxation!
//         if (i!=i_fine) {
//             prolongate_up(grids[i],grids[i+1]);
//         }
//         // reset sor
//         sor = 1.;
//     }
// }






/////////////////////////////
// UP CYCLE 

void Multigrid::up_cycle(const double epsilon) {
    std::cout << std::setprecision(8);

    auto hline = [](char c, int n) { std::cout << std::string(n, c) << "\n"; };

    hline('=', 52);
    std::cout << "  UP CYCLE\n";
    hline('=', 52);

    flush();
    double diff  = 100.;
    int    count = 0;
    double psi4_max_global = 1.;
    double vmax_global = 0.;
    double rhomax_global = 0.;
    double cfl_mod = 1.;

    for (size_t i = i_coarse; i <= i_fine; i++)
    {
        grids[i].fill_all_ghosts();
        count = 0;

        hline('-', 52);
        std::cout << "  GRID " << i << " / " << i_fine;
        if      (i == i_coarse) std::cout << "   [coarse]";
        else if (i == i_fine)   std::cout << "   [fine]";
        else                    std::cout << "   [intermediate]";
        std::cout << "\n";
        hline('-', 52);

        while (true) {

            // multiple consecutive relaxations
            for (int q = 0; q < p.iter; q++) 
            {
                psi4_max_global = std::pow(grids[i].max_psi(),4);
                rhomax_global = grids[i].m_source.max_rho();
                vmax_global = grids[i].m_source.max_v();
                cfl_mod = rhomax_global * 2. * M_PI 
                        * 5. * psi4_max_global;
                cfl_mod = std::max(1.,cfl_mod);
                cfl_mod = 1.;
                grids[i].m_UR = 1./cfl_mod;
                grids[i].relax();
            }

            // double stiffness = 10. * M_PI * std::pow(grids[i].max_psi(), 4) * rhomax_global;
            // double dt_used   = grids[i].m_dt * grids[i].m_UR;
            // std::cout << "  stiffness*dt = " << stiffness * dt_used << "\n";

            diff = std::abs(grids[i].m_delta);
            if (diff < epsilon) {
                std::cout << "  converged after " << count * p.iter 
                          << " iterations  [diff = " << diff << "]\n";
                break;
            }

            // progress line
            // bar width scales with log convergence
            int bar_width = 20;
            double frac = std::log10(diff / epsilon);          // >0 means not converged
            frac = std::max(0., std::min(1., 1. - frac / 6.)); // normalise over ~6 decades
            int filled = static_cast<int>(frac * bar_width);
            std::string bar = std::string(filled, '#') + std::string(bar_width - filled, '.');

            // printing to terminal
            double iW  = grids[i].int_W();
            double iP  = grids[i].int_psi();
            std::cout << "  " << count * p.iter
                    << " it  [" << bar << "]"
                    << "  diff "    << std::scientific << std::setw(11) << diff
                    << std::fixed   << std::setprecision(3)
                    << "  cfl "     << p.CFL / cfl_mod
                    << "  cfl_mod " << std::setprecision(5) << cfl_mod
                    << "  intW "    << iW
                    << "  intp "    << iP
                    << "\n";

            count += 1;
        }

        if (i != i_fine) {
            std::cout << "  prolongating to grid " << i + 1 << " ...\n";
            prolongate_up(grids[i], grids[i + 1]);
        }
    }

    hline('=', 52);
    std::cout << "  UP CYCLE COMPLETE\n";
    hline('=', 52);
    std::cout << "\n";
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
        // std::cout << "Relax grid " << i << "\n"; 
        for (int q = 0; q<p.iter; q++)
        {
            grids[i].relax();
        }
        // std::cout << "Restrict " << i << " to " << i-1 <<"\n";
        restrict_down(grids[i],grids[i-1]);
    }

    // Prolongation 

    // loop up over grids to set initial data on all levels
    for (int i = i_bot; i < i_top; i++)
    {
        // std::cout << "Relax grid " << i << "\n"; 
        for (int q = 0; q<p.iter; q++)
        {
            grids[i].relax();
        }
        // std::cout << "Prolong " << i << " to " << i+1 <<"\n";
        prolongate_up(grids[i],grids[i+1]);
    }

    // std::cout << "Fine level relax\n"; 
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
    double valW = 0;
    double valPsi = 0;
    // fine indices
    int i1 = 0;
    int i2 = 0;
    int j1 = 0;
    int j2 = 0;
    #pragma omp parallel for schedule(static) private(i1,i2,j1,j2,valW,valPsi)
    for (int j = coarse.jmin; j < coarse.jmax; j++) {
        for (int i = coarse.imin; i < coarse.imax; i++) {
            // fine indices
            i1 = i*2 - p.ng;
            i2 = i*2 - p.ng + 1;
            j1 = j*2 - p.ng;
            j2 = j*2 - p.ng + 1;

            valW = fine.W.get_data(i1,j1) + 
                    fine.W.get_data(i2,j1) + 
                    fine.W.get_data(i1,j2) + 
                    fine.W.get_data(i2,j2);
            valPsi = fine.psi.get_data(i1,j1) + 
                    fine.psi.get_data(i2,j1) + 
                    fine.psi.get_data(i1,j2) + 
                    fine.psi.get_data(i2,j2);
            valW *= 0.25;
            valPsi *= 0.25;

            coarse.W.set_data(valW,i,j);
            coarse.psi.set_data(valPsi,i,j);
        }
    }
}


void Multigrid::prolongate_up(const Grid &coarse, Grid &fine) {
    double x = 0.;
    double y = 0.;
    double valW = 0.;
    double valPsi = 0.;
    #pragma omp parallel for schedule(static) private(x,y,valW,valPsi)
    for (int j = fine.jmin; j < fine.jmax; j++) {
        for (int i = fine.imin; i < fine.imax; i++) {

            // ask for bilinear interpolation 
            x = fine.W.get_x(i,j);
            y = fine.W.get_y(i,j);
            valW = coarse.W.bilinear_interp_xy(x,y);
            valPsi = coarse.psi.bilinear_interp_xy(x,y);
            fine.W.set_data(valW,i,j);
            fine.psi.set_data(valPsi,i,j);
        }
    }
}
