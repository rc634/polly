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

    // the two pieces of the ADM mass, each integrated over the box
    // r <= Lx, z <= Ly so the box size dependence can be watched.
    // see Multigrid::check_mass_boxes for where the prefactors come from
    double int_W_term(double Lx, double Ly);   // the extrinsic curvature part
    double int_rho_box(double Lx, double Ly);  // the matter part

    // PDE operator F(u) at one live cell, without tau
    void pde_operator(int i, int j, double& F_W, double& F_psi,
                      double& src_W, double& src_psi);

    // fills xi = F(u) - tau on the live cells
    void compute_residual();

    // parabolic diffusion type relaxation
    void relax();

    // true  : W is ODD about the equator, so net J = 0 and the far field leads
    //         with the l = 2 term  W ~ A sin(th) cos(th) / R^3
    // false : W is EVEN about the equator, the usual l = 1 term
    //         W ~ J sin(th) / R^2
    // public so a run can pick its parity before the first solve
    bool m_antisymmW_BC = false;

    // exterior call routine - selectes which private function to call
    void fill_all_ghosts();
private :
    // interior precise BC functions 
    void fill_all_ghosts_robin();
    void fill_all_ghosts_antisymmW();
    void fill_all_ghosts_ADM();
    void fill_all_ghosts_dirichlet();
public :

    // save data 
    void save_data(const std::string& filename);
    void save_state();

    void hello();
};

