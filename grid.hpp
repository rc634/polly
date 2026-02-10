#pragma once
#include <iostream>
#include "field.hpp"
#include "params.hpp"

class Grid {
public:

    // total number of points
    int num_vars = 1;

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

    // set of physics fields 
    Field W;
    Field psi;

    // place in heirarchy
    int grid_n; 

    // params
    Params p;

    // the effective timestep
    double m_dt;

    // error or delta 
    double m_delta;

    Grid();

    void init(int i);

    void ID_gaussian(double x0, double y0, double sig);

    void ID_zeros();

    void ID_pdisc(double R, double density);

    double field_integral();
    double int_W();
    double int_psi();

    // parabolic diffusion type relaxation
    void relax();

    void fill_all_ghosts();

    // save data 
    void save_data(const std::string& filename);
    void save_state();

    void hello();
};

