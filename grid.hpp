#pragma once
#include <iostream>
#include "field.hpp"
#include "params.hpp"

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

    void ID_gaussian(double x0, double y0, double sig);

    // save data 
    void save_data();

    void hello();
};

