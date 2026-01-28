#include <fstream>
#include <string>
#include "grid.hpp"
#include <iomanip> // for precision output
#include <algorithm> // for std::min


Grid::Grid() {
}

void Grid::init(int i) {
    grid_n = i;
    f1.init(grid_n);
    // courant friendly timestep
    m_dt = p.CFL*std::min(f1.dx,f1.dy)*std::min(f1.dx,f1.dy);
}

void Grid::ID_gaussian(double x0, double y0, double sig) {
    for (int j = f1.jmin; j < f1.jmax; j++) {
        for (int i = f1.imin; i < f1.imax; i++) {
            // r^2 from centre of gaussian
            double rr = pow((x0-f1.get_x(i,j)),2) + pow((y0-f1.get_y(i,j)),2);
            // gaussian
            double val = exp(-(rr)/(2.*sig*sig));
            // set data 
            f1.set_data(val,i,j);
        }
    }
}

double Grid::L2norm() {
    double L2 = 0.;
    for (int k = 0; k < f1.n_flat; k++) {
        L2 += f1.data[k]*f1.data[k];
    }
    return sqrt(L2);
}


void Grid::relax() {
    // One Relaxation step 
    
    // loop over all live cells 
    for (int j = f1.jmin; j < f1.jmax; j++) {
        for (int i = f1.imin; i < f1.imax; i++) {
            // time evolution step

            // coordinate 
            double x = f1.get_x(i,j);
            double y = f1.get_y(i,j);
            double r = sqrt(x*x + y*y);

            // gaussian stats
            double x0 = 0.0;
            double y0 = 0.5;
            double r0 = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));
            double rmax = 0.05;
            double sigx = 0.1;
            double sigy = 0.15;

            // source terms
            double src_f1 = 0;
            if (r0 < rmax) {
                // src_f1 = exp(-(x-x0)*(x-x0)/(2*sigx*sigx))
                //        * exp(-(y-y0)*(y-y0)/(2*sigy*sigy));
                src_f1 = 1.;
            }

            // laplacian is zero WARNING!
            double df1_dt = f1.cartesian_laplacian(i,j) - src_f1;

            // tiemstep f_new = f_old + dt * df/dt
            double new_f1 = f1.get_data(i,j) + df1_dt * m_dt;

            // set new-data stage 
            f1.set_new_data(new_f1,i,j);
        }
    }

    // Save over the old data - 
    // - with the finished new data
    f1.save_new_data();

    // fill all ghosts here ?
}



void Grid::save_data() {
    std::string filename = "data/f.dat";
    std::ofstream file(filename);

    if (!file) {
        std::cerr << "Error: could not open file " << filename << "\n";
        return;
    }

    std::cout << "Saving " << filename << std::endl;

    // 8 significant figures (not 8 decimal places!)
    file << std::scientific << std::setprecision(p.save_precision);

    // auto f_ = f1; // lazy 

    // write everything 
    // including ghost cells for now 
    for (int j = 0; j < f1.nyg; j++) {
        for (int i = 0; i < f1.nxg; i++) {
            file << f1.get_data(i,j) ;
            if (i==f1.nxg-1) {
                break;
            }
            file << ",";
        }
        file << "\n";
    }

    file.close();
}

void Grid::hello() {
    std::cout << "Hello from Grid :  " << grid_n << "\n"; 
    f1.hello();
}
