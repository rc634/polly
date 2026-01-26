#include <fstream>
#include <string>
#include "grid.hpp"


Grid::Grid() {
}

void Grid::init(int i) {
    grid_n = i;
    f1.init(grid_n);
}


void Grid::save_data() {
    std::string filename = "f.dat";
    std::ofstream file(filename);

    if (!file) {
        std::cerr << "Error: could not open file " << filename << "\n";
        return;
    }

    std::cout << "Saving " << filename << std::endl;

    // auto f_ = f1; // lazy 

    // write everything 
    // including ghost cells for now 
    for (int j = 0; j < f1.ny; j++) {
        for (int i = 0; i < f1.nx; i++) {
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
