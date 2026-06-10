// params.hpp
#pragma once

struct Params {
    int num_grids = 6;       // grids [0, num_grids-1]
    int nx = 128;             // grid points in x on coarsest level
    int ny = 128;             // grid points in y on coarsest level
    int ng = 2;              // number of ghost cells per edge
    double tol = 1e-15;       // convergence tolerance
    
    // number of 4th order stencil steps to do near end for fine resolution?
    double iter = 256; // relaxation steps after/before prolongation/restriction

    // lewey friedrich stabiilty thingy 
    // normally 0.25, but the W/x^2 term prefers 0.125
    double CFL = 0.125;

    // physical boundaries, upper lower
    double xL = 0.;
    double xU = 40.;
    double yL = 0.;
    double yU = 40.;

    // data saving
    int save_precision = 16;
};

// notes 

// convergence ..
