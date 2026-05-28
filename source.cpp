#include "source.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <iomanip>

#define _USE_MATH_DEFINES

// -----------------------------------------------------------------------
// Init: analytic — evaluate onto grid and store
// -----------------------------------------------------------------------

void Source::init_analytic(int nxg, int nyg, int ng,
                           double dx, double dy,
                           double xL, double yL) {
    m_nxg = nxg;
    m_nyg = nyg;
    m_v.assign(nxg * nyg, 0.);
    m_rho.assign(nxg * nyg, 0.);

    for (int j = 0; j < nyg; j++) {
        for (int i = 0; i < nxg; i++) {
            // cell-centre coordinates (ghost-inclusive, same convention as Field)
            double x = xL + (i - ng + 0.5) * dx;
            double y = yL + (j - ng + 0.5) * dy;

            // fractional coords wrt ellipsoid
            double xf = x / m_width;
            double yf = y / m_height;
            double rf     = std::sqrt(xf*xf + yf*yf);
            // packet encoding overall compact support
            double packet = 0.;
            double GERRY = sqrt(3.);

            if (xf*xf + yf*yf < 1.) {
                // double env    = 1. + std::cos(rf * M_PI);
                // packet = env * env / 4.; //normalise to packet(0)=1
                packet = exp(GERRY*GERRY*rf*rf/(rf*rf-1.));
            }

            m_v[flat(i,j)]   = m_v_0;  // raw v
            m_rho[flat(i,j)] = m_rho_0 * packet;  // raw rho
        }
    }

    std::cout << " - Source: analytic  rho_0 =" << m_rho_0
              << "  omega_0 =" << m_v_0 << "\n";
}


// -----------------------------------------------------------------------
// Init: load from file
// -----------------------------------------------------------------------

void Source::init_from_file(const std::string& W_path, const std::string& psi_path,
                            int nxg, int nyg, int ng) {
    m_nxg = nxg;
    m_nyg = nyg;
    // ng stored for reference but data is read as full live-cell grid (no ghosts in file)
    int nx = nxg - 2*ng;
    int ny = nyg - 2*ng;

    m_v.assign(nxg * nyg, 0.);
    m_rho.assign(nxg * nyg, 0.);

    auto load = [&](const std::string& path, std::vector<double>& table) {
        std::ifstream f(path);
        if (!f) throw std::runtime_error("Source: cannot open " + path);
        std::string line, token;
        int j = ng;
        while (std::getline(f, line) && j < ng + ny) {
            std::istringstream ss(line);
            int i = ng;
            while (std::getline(ss, token, ',') && i < ng + nx) {
                table[flat(i, j)] = std::stod(token);
                i++;
            }
            j++;
        }
        std::cout << " - Source: loaded " << path << "\n";
    };

    load(W_path,   m_v);
    load(psi_path, m_rho);
}


// -----------------------------------------------------------------------
// Get — simple table lookups, used by multigrid flush
// -----------------------------------------------------------------------

double Source::get_v(int i, int j) const {
    return m_v[flat(i, j)];
}

double Source::get_rho(int i, int j) const {
    return m_rho[flat(i, j)];
}

// -----------------------------------------------------------------------
// Get — raw/eff, interpretation depends on EFFECTIVE_SOURCE flags
// -----------------------------------------------------------------------

double Source::get_v_raw(int i, int j, double x, double psi_ij) const {
    if (EFFECTIVE_V) {
        return m_v[flat(i,j)] / std::pow(psi_ij, 5);  // inverse
    } else {
        return m_v[flat(i,j)];
    }
}

double Source::get_v_eff(int i, int j, double x, double psi_ij) const {
    if (EFFECTIVE_V) {
        return m_v[flat(i,j)];
    } else {
        return std::pow(psi_ij, 5) * m_v[flat(i,j)];
    }
}

double Source::get_rho_raw(int i, int j, double x, double psi_ij) const {
    if (EFFECTIVE_RHO) {
        double v_raw = get_v_raw(i, j, x, psi_ij);
        double denom = 1. - v_raw * v_raw * x * x * std::pow(psi_ij, 4);
        return m_rho[flat(i,j)] * denom / std::pow(psi_ij, 5);
    } else {
        return m_rho[flat(i,j)];
    }
}

double Source::get_rho_eff(int i, int j, double x, double psi_ij) const {
    if (EFFECTIVE_RHO) {
        return m_rho[flat(i,j)];
    } else {
        double v_raw = get_v_raw(i, j, x, psi_ij);
        double rho_raw = m_rho[flat(i,j)];
        double denom = 1. - v_raw * v_raw * x * x * std::pow(psi_ij, 4);
        return std::pow(psi_ij, 5) * rho_raw / denom;
    }
}

// -----------------------------------------------------------------------
// Get — full PDE source terms, always uses effective quantities
// -----------------------------------------------------------------------

double Source::get_src_W(int i, int j, double x, double psi_ij) const {
    return 8. * M_PI * x * get_v_eff(i,j,x,psi_ij) * get_rho_eff(i,j,x,psi_ij);
}

double Source::get_src_psi(int i, int j, double x, double psi_ij) const {
    return -2. * M_PI * get_rho_eff(i,j,x,psi_ij);
}

double Source::max_rho() const {
    double rho_max = 0.;
    for (int j = 0; j < m_nyg; j++) {
        for (int i = 0; i < m_nxg; i++) {
            double rho_ij = std::abs(m_rho[flat(i,j)]);
            if (rho_ij > rho_max) rho_max = rho_ij;
        }
    }
    return rho_max;
}

double Source::max_v() const {
    double v_max = 0.;
    for (int j = 0; j < m_nyg; j++) {
        for (int i = 0; i < m_nxg; i++) {
            double v_ij = std::abs(m_v[flat(i,j)]);
            if (v_ij > v_max) v_max = v_ij;
        }
    }
    return v_max;
}


// -----------------------------------------------------------------------
// Save2 — same format as Grid::save_state()
// -

void Source::save(const std::string& out_dir, int ng,
                  const Field& psi_field) const {
    int nx = m_nxg - 2*ng;
    int ny = m_nyg - 2*ng;

    std::ofstream f_rho_raw(out_dir + "/rho-raw.dat");
    std::ofstream f_rho_eff(out_dir + "/rho-eff.dat");
    std::ofstream f_src_psi(out_dir + "/source-psi.dat");
    std::ofstream f_v_raw  (out_dir + "/v-raw.dat");
    std::ofstream f_v_eff  (out_dir + "/v-eff.dat");
    std::ofstream f_src_W  (out_dir + "/source-W.dat");

    // set precision on all
    for (auto* f : {&f_rho_raw, &f_rho_eff, &f_src_psi,
                    &f_v_raw,   &f_v_eff,   &f_src_W})
        *f << std::scientific << std::setprecision(10);

    for (int j = ng; j < ng + ny; j++) {
        for (int i = ng; i < ng + nx; i++) {
            double x      = psi_field.get_x(i, j);
            double psi_ij = psi_field.get_data(i, j);

            std::string sep = (i < ng + nx - 1) ? "," : "\n";

            f_rho_raw << get_rho_raw(i,j,x,psi_ij) << sep;
            f_rho_eff << get_rho_eff(i,j,x,psi_ij) << sep;
            f_src_psi << get_src_psi(i,j,x,psi_ij) << sep;
            f_v_raw   << get_v_raw(i,j,x,psi_ij)   << sep;
            f_v_eff   << get_v_eff(i,j,x,psi_ij)   << sep;
            f_src_W   << get_src_W(i,j,x,psi_ij)   << sep;


        }
    }

    std::cout << "Saved 6 source files to " << out_dir << "\n";
}

