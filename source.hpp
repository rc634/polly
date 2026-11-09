#pragma once
#include <string>
#include <vector>
#include "field.hpp" // used for saving sources mixed with geometry

class Source {
public:
    Source() {}

    static constexpr bool ANALYTIC = true;
    static constexpr bool EFFECTIVE_RHO = true;
    static constexpr bool EFFECTIVE_V = true;

    // smooth bump profiles
    void init_bump(int nxg, int nyg, int ng,
                   double dx, double dy,
                   double xL, double yL);

    // flat density profile nakamura like profiles 
    void init_naka(int nxg, int nyg, int ng,
                   double dx, double dy,
                   double xL, double yL);

    // symmetric jacobi profiles 
    void init_SJ(int nxg, int nyg, int ng,
                 double dx, double dy,
                 double xL, double yL);

    // Init by loading source-W.dat and source-psi.dat
    void init_from_file(const std::string& W_path, const std::string& psi_path,
                        int nxg, int nyg, int ng);


    // stored quantities (interpretation depends on EFFECTIVE_SOURCE flag)
    double get_v(int i, int j) const;
    double get_rho(int i, int j) const;

    // raw physical quantities
    double get_v_raw(int i, int j, double x, double psi_ij) const;
    double get_rho_raw(int i, int j, double x, double psi_ij) const;

    // effective quantities
    double get_v_eff(int i, int j, double x, double psi_ij) const;
    double get_rho_eff(int i, int j, double x, double psi_ij) const;

    // full PDE source terms
    double get_src_W(int i, int j, double x, double psi_ij) const;
    double get_src_psi(int i, int j, double x, double psi_ij) const;

    // setter functions
    void set_v(int i, int j, double val)   { m_v[flat(i,j)]   = val; }
    void set_rho(int i, int j, double val) { m_rho[flat(i,j)] = val; }

    // ADM mass from integrating effective source over flat space
    // M_ADM = 2π ∫∫ ρ_eff · x dx dy   (exact for conformally flat Hamiltonian constraint)
    double m_M_ADM = 0.;

    // max value finders for courant
    double max_rho() const;
    double max_v() const;

    // Save to dat
    void save(const std::string& out_dir, int ng, const Field& psi_field) const;

    // post-solve check that the matter is still physical. the raw density
    // carries a factor 1 - v_phys^2 with v_phys = v_raw x psi^2, which turns
    // negative once rotation drives the matter superluminal. needs the solved
    // psi, so it can only be run after the multigrid solve.
    // returns true if everything stayed positive
    bool sanity_check(int ng, const Field& psi_field) const;

    void set_source_params(double rho_0, double v_0, double height, double width) {
        m_rho_0  = rho_0;
        m_v_0    = v_0;
        m_height = height;
        m_width  = width;
    }

private:

// arrays 
    std::vector<double> m_v; // source rotation array
    std::vector<double> m_rho; // source density array

    // matter density stuff 
    double m_rho_0  =  12.0;// 4.0; // 3.4;
    double m_v_0    = 0.4;
    // bump steepness, exp(sigma^2 rf^2 / (rf^2 - 1)) in init_bump
    double m_sigma   = std::sqrt(3.);

    // matter shape 
    double m_height = 4.; // 0.5;
    double m_width  = 0.5; // sqrt(2.);

    // rotation stuff
    // v = (v_0/(x_f + eps_v)) (1 + a_v * bump(sigma_v)). a_v = 0 is uniform Omega,
    // positive spins the core up relative to the crust
    double m_sigma_v = std::sqrt(2.);
    double m_a_v = 3.0;
    // regulator on the 1/x_f profile, so the linear velocity Omega x goes to zero
    // on the axis instead of to a constant. in units of the width, and far below
    // the finest dx/2, so the gridded data only moves at the 1e-8 level
    double m_eps_v = 1e-10;

    // initial data binary separation ratio, make > 1
    double id_hsep_r = 2.0;

    // grid stuff
    int m_nxg = 0;
    int m_nyg = 0;
    int m_ng  = 0;
    double m_dx = 0., m_dy = 0., m_xL = 0.;
    int flat(int i, int j) const { return j * m_nxg + i; }

    // compact support bump, 1 at the centre and 0 outside rf = 1. shared by
    // the density packet and the rotation profile, each with its own sigma
    static double bump(double rf, double sigma) {
        if (rf >= 1.) return 0.;
        return std::exp(sigma*sigma*rf*rf/(rf*rf - 1.));
    }

    static double heaviside(double rf) {
        if (rf >= 1.) return 0.;
        return 1.;
    }


    void integrate_M_ADM();
};



