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

    // Init by evaluating analytic source onto grid
    // Requires ghost-inclusive grid dimensions and cell-centre coordinate info
    void init_analytic(int nxg, int nyg, int ng,
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

    // max value finders for courant
    double max_rho() const;
    double max_v() const;

    // Save to dat
    void save(const std::string& out_dir, int ng, const Field& psi_field) const;

    void set_source_params(double rho_0, double v_0, double height, double width) {
        m_rho_0  = rho_0;
        m_v_0    = v_0;
        m_height = height;
        m_width  = width;
    }

private:
    double m_rho_0  = 4.0;
    double m_v_0    = 1.0;
    double m_height = 0.5;
    double m_width  = sqrt(2.);
    std::vector<double> m_v; // source rotation array
    std::vector<double> m_rho; // source density array
    int m_nxg = 0;
    int m_nyg = 0;
    int flat(int i, int j) const { return j * m_nxg + i; }
};



