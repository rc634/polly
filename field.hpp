#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include "params.hpp"

class Field {
public:
    Params p;

    // grid number (heirarchy)
    int grid_n; // 0 minimum

    // physical problem size 
    int nx;
    int ny;

    // numerical problem size 
    int nxg;
    int nyg;

    // ghosts 
    int ng = p.ng;

    std::vector<double> data; // flattened 2D array

    Field();

    void init(int i);

    void set_data(double val, int i, int j);

    double get_data(int i, int j);

    // index getter
    inline int index(int i, int j);

    void hello();
};