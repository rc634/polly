// params.hpp
#pragma once

struct Params {
    int num_grids = 3;       // number of levels in the multigrid hierarchy
    int nx = 16;             // grid points in x on coarsest level
    int ny = 16;             // grid points in y on coarsest level
    int ng = 2;              // number of ghost cells per edge
    double tol = 1e-6;       // convergence tolerance

    // physical boundaries
    double xL = 0.;
    double xR = 1.;
    double yL = 0.;
    double yR = 1.;
};
