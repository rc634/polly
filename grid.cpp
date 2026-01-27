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
    m_dt = p.CFL*std::min(f1.dx,f1.dy);
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



void Grid::relax() {
    // One Relaxation step 
    
    // loop over all live cells 
    for (int j = f1.jmin; j < f1.jmax; j++) {
        for (int i = f1.imin; i < f1.imax; i++) {
            // time evolution step

            // laplacian is zero WARNING!
            double df1_dt = f1.cartesian_laplacian(i,j);

            // tiemstep f_new = dt * df/dt + f
            double new_f1 = df1_dt * m_dt + f1.get_data(i,j);

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
