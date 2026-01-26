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

    // constructor
    Multigrid();

    void init();

    // member functions
    void vcycle();

    // save data 
    void save_data();

    // hello world
    void hello();
};
