// params.hpp
#pragma once

struct Params {
    int num_grids = 4;       // number of levels in the multigrid hierarchy, min 1
    int nx = 32;             // grid points in x on coarsest level
    int ny = 32;             // grid points in y on coarsest level
    int ng = 2;              // number of ghost cells per edge
    double tol = 1e-6;       // convergence tolerance
    
    // number of 4th order stencil steps to do near end for fine resolution?
    double smooth_steps = 100; 

    // physical boundaries, upper lower
    double xL = 0.;
    double xU = 1.;
    double yL = 0.;
    double yU = 1.;

    // data saving
    int save_precision = 3;
};
