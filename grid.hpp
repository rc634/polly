#pragma once
#include <iostream>
#include "field.hpp"
#include "params.hpp"
#include "source.hpp"

class Grid {
public:

    // total number of points
    int num_vars = 1;

    // gridpoints
    int nx;
    int ny;
    double dx;
    double dy;

    // physical sizes including ghosts
    double xL, xU, yL, yU;

    // numerical problem size 
    int nxg;
    int nyg;
    int n_flat;

    // used for looping physics 
    // for (i = imin; i < imax; i++ )...
    int imin;
    int imax;
    int jmin;
    int jmax;

    // ghost layers (e.g. per side)
    int ng;

    // set of physics fields 
    Field W;
    Field psi;

    // matter source object
    Source m_source;

    // place in heirarchy
    int grid_n; 

    // params
    Params p;

    // the effective timestep
    double m_dt;

    // currently unused under-relaxation
    // under-relaxation, to account for sources
    double m_UR;

    // volume weighted RMS of the elliptic residual, and of the source in the
    // same norm. both are set by every relax()
    double m_res = 0.;
    double m_src_norm = 0.;

    // relative residual, the dimensionless quantity convergence tests on.
    // with v_0 = 0 the W source vanishes, but psi's does not, so the combined
    // denominator stays alive. falls back to the absolute residual for the
    // degenerate case of no source at all
    double rel_res() const {
        return (m_src_norm > 1e-300) ? m_res / m_src_norm : m_res;
    }

    Grid();

    void init(int i);

    void ID_gaussian(double x0, double y0, double sig);

    void ID_CappedSchwarzschild();

    void ID_pdisc(double R, double density);

    // returns max value of psi for courant time step stuff
    double max_psi();

    double field_integral();
    double int_W();
    double int_psim1();

    // parabolic diffusion type relaxation
    void relax();

    // exterior call routine - selectes which private function to call
    void fill_all_ghosts();
private :
    // interior precise BC functions 
    void fill_all_ghosts_robin();
    void fill_all_ghosts_ADM();
    void fill_all_ghosts_dirichlet();
public :

    // save data 
    void save_data(const std::string& filename);
    void save_state();

    void hello();
};

