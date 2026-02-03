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

    std::vector<double> data; // flattened 2D array
    std::vector<double> data_new; // flattened 2D array

    Field();

    void hello();

    void init(int i);

    void set_data(double val, int i, int j);

    void set_new_data(double val, int i, int j);

    void save_new_data();

    void copy_old_ghosts_to_new();

    // getters 

    double get_new_data(int i, int j);

    double get_data(int i, int j) const;

    double get_x(int i, int j) const;

    double get_y(int i, int j) const;

    inline int index(int i, int j) const;

    // bilinear interpolator

    double bilinear_interp_xy(double x, double y) const;

    // L2 norm new - old data

    double delta_data();

    // differential operators

    double d1x(int i, int j);
    double d1y(int i, int j);
    double d2x(int i, int j);
    double d2y(int i, int j);
    double cartesian_laplacian(int i, int j);
    double cylindrical_laplacian(int i, int j, double x);
    double cylindrical_laplacian_bad(int i, int j, double x);
};