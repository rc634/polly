#include "field.hpp"
#include "stencils.hpp"


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
    // flattened array index
    n_flat = nxg*nyg;

    // data arrays
    data.resize(n_flat,0.0);
    data_new.resize(n_flat,0.0);
    rhs.resize(n_flat,0.0);
    tau.resize(n_flat,0.0);

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

void Field::set_new_data(double val, int i, int j) {
    data_new[i + nxg * j] = val;
}

double Field::get_new_data(int i, int j) {
    return data_new[i + nxg * j];
}

double Field::delta_data() {
    // loop all cells including ghosts
    // flattened array type loop
    double out = 0.0;
    for (int k = 0; k < n_flat; k++) {
        out += pow(data[k] - data_new[k],2);
    }
    return out;
}

void Field::save_new_data() {
    // loop all cells including ghosts
    // flattened array type loop
    for (int k = 0; k < n_flat; k++) {
        data[k] = data_new[k];
    }
}

// copy ghosts from old data to new data 
void Field::copy_old_ghosts_to_new() {
    // two ghosts thickness hard coded!

    // loop left boundary
    for (int j = 0; j < nyg; j++) {
        set_new_data(get_data(0,j),0,j);
        set_new_data(get_data(1,j),1,j);
    }
  
    // loop right boundary
    for (int j = 0; j < nyg; j++) {
        set_new_data(get_data(nxg-2,j), nxg-2, j);
        set_new_data(get_data(nxg-1,j), nxg-1, j);
    }
    // loop top boundary
    for (int i = 0; i < nxg; i++) {
        set_new_data(get_data(i,nyg-2), i, nyg-2);
        set_new_data(get_data(i,nyg-1), i, nyg-1);
    }
    // loop bottom boundary
    for (int i = 0; i < nxg; i++) {
        set_new_data(get_data(i,1), i,1);
        set_new_data(get_data(i,0), i,0);
    }
}

void Field::hello() {
    std::cout << "Hello from Field : " << nxg << " x " << nyg << "\n";
}

// calculus 
double Field::d1x(int i, int j) {
    // stencil loading 
    const auto& stencil = fd::d1_central_2;

    // perform stencil deriv
    double deriv = 0.0;
    for (std::size_t k = 0; k < stencil.offset.size(); k++) {
        // create flattened index 
        int m = index(i + stencil.offset[k],j);
        deriv += stencil.coeff[k] * data[m];
    }
    return deriv / dx;
}

double Field::d1y(int i, int j) {
    // stencil loading 
    const auto& stencil = fd::d1_central_2;

    // perform stencil deriv
    double deriv = 0.0;
    for (std::size_t k = 0; k < stencil.offset.size(); k++) {
        // create flattened index 
        int m = index(i,j + stencil.offset[k]);
        deriv += stencil.coeff[k] * data[m];
    }
    return deriv / dy;
}

double Field::d2x(int i, int j) {
    // stencil loading 
    const auto& stencil = fd::d2_central_2;

    // perform stencil deriv
    double deriv = 0.0;
    for (std::size_t k = 0; k < stencil.offset.size(); k++) {
        // create flattened index 
        int m = index(i + stencil.offset[k],j);
        deriv += stencil.coeff[k] * data[m];
    }
    return deriv / (dx * dx);
}

double Field::d2y(int i, int j) {
    // stencil loading 
    const auto& stencil = fd::d2_central_2;

    // perform stencil deriv
    double deriv = 0.0;
    for (std::size_t k = 0; k < stencil.offset.size(); k++) {
        // create flattened index 
        int m = index(i,j + stencil.offset[k]);
        deriv += stencil.coeff[k] * data[m];
    }
    return deriv / (dy * dy);
}

double Field::cartesian_laplacian(int i, int j) {
    return d2x(i,j) + d2y(i,j);
}

double Field::cylindrical_laplacian(int i, int j, double x) {
    double lap_y = d2y(i,j); // cartesian like z coord
    double f0 = get_data(i-1,j);
    double f1 = get_data(i,j);
    double f2 = get_data(i+1,j);
    double lap_x_curvilinear = ((x+0.5*dx)*(f2-f1) - (x-0.5*dx)*(f1-f0))/(x*dx*dx);
    return lap_y + lap_x_curvilinear;
}

double Field::cylindrical_laplacian_bad(int i, int j, double x) {
    return d2x(i,j) + d2y(i,j) + d1x(i,j)/x;
}