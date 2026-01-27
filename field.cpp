#include "field.hpp"


Field::Field() {
}

void Field::init(int i) {
    // grid specific
    grid_n = i;
    nx = (int) pow(2,grid_n) * p.nx;
    ny = (int) pow(2,grid_n) * p.ny;
    // cell centered
    dx = (p.xU-p.xL)/nx;
    dy = (p.yU-p.yL)/ny;

    // deal with ghosts 
    ng = p.ng;
    nxg = nx + 2*ng;
    nyg = ny + 2*ng;
    data.resize(nxg*nyg,0.0);

    // indices 
    imin = ng;
    imax = ng + nx;
    jmin = ng;
    jmax = ng + ny;

}

// index getter
inline int Field::index(int i, int j) const {
    return i + nxg * j;
}

// data getter
double Field::get_data(int i, int j) const {
    return data[i + nxg * j];
}

// also do get x get y with the single flattened index?
// write index UN flatteners? e.g. k -> i,j

// cell centred coordinates!
double Field::get_x(int i, int j) const {
    return (i - p.ng + 0.5)*dx + p.xL;
}

double Field::get_y(int i, int j) const {
    return (j - p.ng + 0.5)*dy + p.yL;
}

double Field::bilinear_interp_xy(double x, double y) const {
    double val = 0.;
    double ireal = (x-p.xL)/dx + p.ng - 0.5; 
    double jreal = (y-p.yL)/dy + p.ng - 0.5; 
    int iL = (int) floor(ireal);
    int jL = (int) floor(jreal);
    double a = ireal - iL;
    double b = jreal - jL;

    // gather data from 4 nearest neighbours 
    double f00 = get_data(iL,jL); 
    double f01 = get_data(iL,jL+1); 
    double f10 = get_data(iL+1,jL); 
    double f11 = get_data(iL+1,jL+1);

    // bilinear interpolation
    val = (1-a)*(1-b)*f00 
        + a*(1-b)*f10
        + b*(1-a)*f01
        + a*b*f11;

    return val;
}

void Field::set_data(double val, int i, int j) {
    data[i + nxg * j] = val;
}

void Field::hello() {
    std::cout << "Hello from Field : " << nxg << " x " << nyg << "\n";
}