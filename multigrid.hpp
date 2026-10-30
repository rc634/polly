#pragma once
#include <iostream>
#include <vector>
#include "grid.hpp"
#include "params.hpp"

class Multigrid {
public:
    Params p;

    // set of grids
    std::vector<Grid> grids;

    // pointer to finest grid
    Grid* fine_grid_ptr;
    Grid* coarse_grid_ptr;

    // info about grids 
    int i_fine, i_coarse; // integers to index fine and coarse grids

    // constructor
    Multigrid();

    void init();

    // save data 
    void save_data(const std::string& filename);

    // set field data pre-solve
    void initial_data();

    // set field data pre-solve
    void analytic_solution();

    // fill ghosts 
    void fill_all_ghosts();

    // flush data down from fine level to coarse and back up to fine again
    void flush();

    // flush data down from fine level to coarse and back up to fine again
    void v_cycle(int i_top, int i_bot);

    // flush data down from fine level to coarse and back up to fine again
    void w_cycle();

    // starts at low res then solves its way up with the threshold
    void up_cycle(const double epsilon);

    // repeated v cycles on the full heirarchy until the fine grid meets the threshold
    void v_solve(const double epsilon);

    // single relax step of fine grid
    void refine();

    // hello world
    void hello();

    // restriction : fine -> coarse
    // promises not to modify fine 
    void restrict_down(const Grid &fine, Grid &coarse);

    // prolongate : coarse -> fine
    // promises not to modify coarse
    void prolongate_up(const Grid &coarse, Grid &fine);


    // FAS restriction : fine -> coarse
    // coarse data = R u_f, and coarse tau = F_c(R u_f) - R xi_f
    void restrict_fas(Grid &fine, Grid &coarse);

    // FAS prolongation : coarse -> fine
    // u_f += P(u_c - R u_f). leaves R u_f in coarse data
    void prolongate_fas(Grid &coarse, Grid &fine);

    // zero tau on every level, so plain relaxation solves the original PDE
    void clear_tau();

    // sources
    void restrict_source_down(const Grid& fine, Grid& coarse);

    // is the converged matter still physical? run on the finest grid at the end
    // of every solve, since it needs psi. returns true if all is well
    bool check_source();

    // the ADM mass split into its matter and extrinsic curvature parts, each
    // integrated over nested boxes, so the box size dependence is visible.
    // returns the full box total, and leaves the three pieces below for main
    // to write into params.dat
    double check_mass_boxes(int n_boxes = 4);

    // full box masses, set by check_mass_boxes at the end of every solve.
    // m_M_ADM is the one downstream tools want: the matter term alone is short
    // by the whole W contribution, which the spinning runs care about
    double m_M_ADM = 0., m_M_rho = 0., m_M_W = 0.;
};
