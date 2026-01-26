#include "field.hpp"


Field::Field() {
}

void Field::init(int i) {
    grid_n = i;
    nx = (int) pow(2,grid_n) * p.nx;
    ny = (int) pow(2,grid_n) * p.ny;
    nxg = nx + 2*ng;
    nyg = ny + 2*ng;
    data.resize(nxg*nyg,0.0);
}

// index getter
inline int Field::index(int i, int j) {
    return i + nxg * j;
}

// data getter
double Field::get_data(int i, int j) {
    return data[i + nxg * j];
}

void Field::set_data(double val, int i, int j) {
    data[i + nxg * j] = val;
}

void Field::hello() {
    std::cout << "Hello from Field : " << nxg << " x " << nyg << "\n";
}