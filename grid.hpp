#pragma once
#include <iostream>
#include "field.hpp"

class Grid {
public:

    // total number of points
    int num_vars = 1;

    // set of physics fields 
    Field f1;

    // place in heirarchy
    int grid_n; 

    Grid();

    void init(int i);

    // save data 
    void save_data();

    void hello();
};

