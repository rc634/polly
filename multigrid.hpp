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

    // constructor
    Multigrid();

    void init();

    // save data 
    void save_data();

    // set field data pre-solve
    void initial_data();

    // fill ghosts 
    void fill_all_ghosts();

    // flush data down from fine level to coarse and back up to fine again
    void flush();

    // flush data down from fine level to coarse and back up to fine again
    void v_cycle();

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
};
