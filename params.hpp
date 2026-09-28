// params.hpp
#pragma once
#include <string>

struct Params {
    int num_grids = 3;       // grids [0, num_grids-1]
    int nx = 256; // 128; // 512;             // grid points in x on coarsest level
    int ny = 256; // 128;// 512;             // grid points in y on coarsest level
    int ng = 2;              // number of ghost cells per edge

    // convergence tolerance on the RELATIVE RESIDUAL ||r|| / ||src||
    // independent of both resolution and dt
    double res_tol = 1e-5;   // 10e-8

    // number of 4th order stencil steps to do near end for fine resolution?
    double iter = 256; // relaxation steps after/before prolongation/restriction

    // lewey friedrich stabiilty thingy 
    // normally 0.25, but the W/x^2 term prefers 0.125
    double CFL = 0.125;

    // physical boundaries, upper lower
    double xL = 0.;
    double xU = 80.; // 160.;
    double yL = 0.;
    double yU = 80.; // 160.;

    // data saving
    int save_precision = 16;
    std::string out_dir = "data";
};

// notes 

// convergence ..
