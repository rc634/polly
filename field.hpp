#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include "params.hpp"

class Field {
public:
    Params p;

    // grid number (heirarchy)
    int grid_n; // 0 minimum

    // gridpoints
    int nx;
    int ny;
    double dx;
    double dy;

    // physical sizes including ghosts
    double xL, xU, yL, yU;

    // numerical problem size 
    int nxg;
    int nyg;
    int n_flat;

    // used for looping physics 
    // for (i = imin; i < imax; i++ )...
    int imin;
    int imax;
    int jmin;
    int jmax;

    // ghost layers (e.g. per side)
    int ng;

    // flattened 2D data arrays (BIG)
    // dynamic solution
    std::vector<double> data; 

    // delete later 
    std::vector<double> data_new; 

    // residual of the PDE, F(u) - tau, evaluated at the current data
    // restricted down to build the coarse grid tau
    std::vector<double> xi;

    // FAS source offset, fixed for a whole v cycle
    // zero on the finest grid, set on restriction to the coarse grid:
    //   tau_c = F_c(R u_f) - R xi_f
    // relaxation drives F(u) - tau -> 0
    std::vector<double> tau;

    Field();

    void hello();

    void init(int i);

    void set_data(double val, int i, int j);

    void set_new_data(double val, int i, int j);

    void save_new_data();

    void copy_old_ghosts_to_new();

    void set_xi(double val, int i, int j);

    void set_tau(double val, int i, int j);

    // getters

    double get_xi(int i, int j) const;

    double get_tau(int i, int j) const;

    double get_new_data(int i, int j);

    double get_data(int i, int j) const;

    double get_x(int i, int j) const;

    double get_y(int i, int j) const;

    inline int index(int i, int j) const;

    // bilinear interpolator

    double bilinear_interp_xy(double x, double y) const;

    // differential operators

    double d1x(int i, int j);
    double d1y(int i, int j);
    double d2x(int i, int j);
    double d2y(int i, int j);
    double cartesian_laplacian(int i, int j);
    double cylindrical_laplacian(int i, int j, double x);
    double cylindrical_laplacian_HC(int i, int j, double x);
};