#include <fstream>
#include <string>
#include "grid.hpp"
#include <iomanip> // for precision output


Grid::Grid() {
}

void Grid::init(int i) {
    grid_n = i;
    f1.init(grid_n);
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



void Grid::save_data() {
    std::string filename = "data/f.dat";
    std::ofstream file(filename);

    if (!file) {
        std::cerr << "Error: could not open file " << filename << "\n";
        return;
    }

    std::cout << "Saving " << filename << std::endl;

    // 8 significant figures (not 8 decimal places!)
    file << std::scientific << std::setprecision(1);

    // auto f_ = f1; // lazy 

    // write everything 
    // including ghost cells for now 
    for (int j = 0; j < f1.nyg; j++) {
        for (int i = 0; i < f1.nxg; i++) {
            file << f1.get_data(i,j) << ",";
        }
        file << "\n";
    }

    file.close();
}

void Grid::hello() {
    std::cout << "Hello from Grid :  " << grid_n << "\n"; 
    f1.hello();
}
