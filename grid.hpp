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

    // params
    Params p;

    // the effective timestep
    double m_dt;

    // error or delta 
    double m_delta;

    Grid();

    void init(int i);

    void ID_gaussian(double x0, double y0, double sig);

    void ID_zeros();

    void ID_pdisc(double R, double density);

    double L2norm();

    double field_integral();

    // parabolic diffusion type relaxation
    void relax();

    void fill_all_ghosts();

    // save data 
    void save_data(const std::string& filename);

    void hello();
};

