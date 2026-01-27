// params.hpp
#pragma once

struct Params {
    int num_grids = 3;       // number of levels in the multigrid hierarchy, min 1
    int nx = 8;             // grid points in x on coarsest level
    int ny = 8;             // grid points in y on coarsest level
    int ng = 2;              // number of ghost cells per edge
    double tol = 1e-6;       // convergence tolerance

    // physical boundaries, upper lower
    double xL = 0.;
    double xU = 1.;
    double yL = 0.;
    double yU = 1.;
};
