// params.hpp
#pragma once

struct Params {
    int num_grids = 5;       // grids [0, num_grids-1]
    int nx = 64;             // grid points in x on coarsest level
    int ny = 64;             // grid points in y on coarsest level
    int ng = 2;              // number of ghost cells per edge
    double tol = 1e-16;       // convergence tolerance
    
    // number of 4th order stencil steps to do near end for fine resolution?
    double iter = 1024; // relaxation steps after/before prolongation/restriction

    // lewey friedrich stabiilty thingy 
    // normally 0.25, but the W/x^2 term prefers 0.125
    double CFL = 0.125;

    // physical boundaries, upper lower
    double xL = 0.;
    double xU = 20.;
    double yL = 0.;
    double yU = 20.;

    // data saving
    int save_precision = 16;
};

// notes 

// convergence ..
