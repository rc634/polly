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


    // sources
    void restrict_source_down(const Grid& fine, Grid& coarse);
};
