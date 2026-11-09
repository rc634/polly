#include "source.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <algorithm>
#include <locale>

#define _USE_MATH_DEFINES

// -----------------------------------------------------------------------
// Init: analytic — evaluate onto grid and store
// -----------------------------------------------------------------------

void Source::init_bump(int nxg, int nyg, int ng,
                       double dx, double dy,
                       double xL, double yL) {
    m_nxg = nxg;
    m_nyg = nyg;
    m_ng  = ng;
    m_dx  = dx;
    m_dy  = dy;
    m_xL  = xL;
    m_v.assign(nxg * nyg, 0.);
    m_rho.assign(nxg * nyg, 0.);

    for (int j = 0; j < nyg; j++) {
        for (int i = 0; i < nxg; i++) {
            double x = xL + (i - ng + 0.5) * dx;
            double y = yL + (j - ng + 0.5) * dy;

            double xf = x / m_width;
            double yf = y / m_height  - id_hsep_r;
            // spheroidal radius
            double rf = std::sqrt(xf*xf + yf*yf);

            // rigid rotation
            // m_v[flat(i,j)]   = m_v_0;

            // differential velocity
            // Omega ~ 1/(x_f + eps_v) rather than 1/x_f, so the linear velocity
            // Omega x goes to zero on the axis instead of to a constant. eps_v
            // is far below the first cell at dx/2, so the gridded data only
            // moves at the 1e-8 level (J_ADM by 3e-10 on paper-diffrotsph-HV/v12)
            m_v[flat(i,j)]   = m_v_0 * (1. + m_a_v * bump(rf, m_sigma_v)) / (xf + m_eps_v);
            
            // matter density // DONT UNCOMMENT!
            m_rho[flat(i,j)] = m_rho_0 * bump(rf, m_sigma);
        }
    }

    integrate_M_ADM();
    std::cout << " - Source: bump  rho_0=" << m_rho_0
              << "  omega_0=" << m_v_0 << "  M_ADM=" << m_M_ADM << "\n";
}


void Source::init_naka(int nxg, int nyg, int ng,
                       double dx, double dy,
                       double xL, double yL) {
    m_nxg = nxg;
    m_nyg = nyg;
    m_ng  = ng;
    m_dx  = dx;
    m_dy  = dy;
    m_xL  = xL;
    m_v.assign(nxg * nyg, 0.);
    m_rho.assign(nxg * nyg, 0.);

    for (int j = 0; j < nyg; j++) {
        for (int i = 0; i < nxg; i++) {
            double x = xL + (i - ng + 0.5) * dx;
            double y = yL + (j - ng + 0.5) * dy;

            double xf = x / m_width;
            double yf = y / m_height;
            // spheroidal radius
            double rf = std::sqrt(xf*xf + yf*yf);

            // vortex velocity
            // same eps_v regulator as init_bump, so Omega x vanishes on the axis
            m_v[flat(i,j)]   = m_v_0 * (1. + m_a_v * heaviside(rf)) / (xf + m_eps_v);
            m_rho[flat(i,j)] = m_rho_0 * heaviside(rf);
        }
    }

    integrate_M_ADM();
    std::cout << " - Source: bump  rho_0=" << m_rho_0
              << "  omega_0=" << m_v_0 << "  M_ADM=" << m_M_ADM << "\n";
}


void Source::init_SJ(int nxg, int nyg, int ng,
                     double dx, double dy,
                     double xL, double yL) {
    m_nxg = nxg;
    m_nyg = nyg;
    m_ng  = ng;
    m_dx  = dx;
    m_dy  = dy;
    m_xL  = xL;
    m_v.assign(nxg * nyg, 0.);
    m_rho.assign(nxg * nyg, 0.);

    for (int j = 0; j < nyg; j++) {
        for (int i = 0; i < nxg; i++) {
            double x = xL + (i - ng + 0.5) * dx;
            double y = yL + (j - ng + 0.5) * dy;

            double xf = x / m_width;
            double yf = y / m_height;
            double rf2 = xf*xf + yf*yf;
            double packet = 0.;

            if (rf2 < 1.) {
                double t = 1. - rf2;
                packet = t*t*t*t;  // (1 - rf^2)^4
            }

            m_v[flat(i,j)]   = m_v_0 * (1. + m_a_v * bump(xf, m_sigma_v));
            m_rho[flat(i,j)] = m_rho_0 * packet;
        }
    }

    integrate_M_ADM();
    std::cout << " - Source: SJ  rho_0=" << m_rho_0
              << "  omega_0=" << m_v_0 << "  M_ADM=" << m_M_ADM << "\n";
}


// -----------------------------------------------------------------------
// M_ADM: exact for conformally-flat Hamiltonian constraint
//   ∇²ψ = -2π ρ_eff   =>   M = 2π ∫∫ ρ_eff x dx dy   (Gauss, axisymmetry)
// -----------------------------------------------------------------------
void Source::integrate_M_ADM() {
    double M_half = 0.;
    for (int j = m_ng; j < m_nyg - m_ng; j++) {
        for (int i = m_ng; i < m_nxg - m_ng; i++) {
            double x = m_xL + (i - m_ng + 0.5) * m_dx;
            M_half += 2. * M_PI * get_rho_eff(i, j, x, 1.) * x * m_dx * m_dy;
        }
    }
    // doubled as integrator uses summetry plane cutting solution in half 
    m_M_ADM = 2 * M_half;
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
// Sanity check — is the converged solution still physical?
//
// the rotating dust carries a factor
//   denom = 1 - v_raw^2 x^2 psi^4 = 1 - v_phys^2
// where v_phys = v_raw x psi^2 is the proper rotational velocity. it turns
// negative as soon as the matter is driven superluminal, and takes the physical
// density with it: with EFFECTIVE_RHO on it is rho_raw that flips, with it off
// rho_eff does. either way the test is the same single condition, v_phys > 1,
// so that is all this reports.
//
// the solver never sees any of it. the elliptic source it relaxes against stays
// positive by construction, so the multigrid converges perfectly happily on data
// that is already unphysical. psi only exists at the end, so the check only
// means anything here.
//
// only cells that actually hold matter are counted. the vacuum carries the same
// factor but multiplies it by zero density, so counting it would report the
// maximum from somewhere the matter never reaches
// -----------------------------------------------------------------------
bool Source::sanity_check(int ng, const Field& psi_field) const {
    int nx = m_nxg - 2*ng;
    int ny = m_nyg - 2*ng;

    double max_vphys = 0.;
    double bad_x = 0., bad_y = 0.;
    long   n_matter = 0, n_bad = 0;

    for (int j = ng; j < ng + ny; j++) {
        for (int i = ng; i < ng + nx; i++) {
            // the bump has compact support, so this is the a,b matter ellipse
            if (m_rho[flat(i,j)] <= 0.) continue;
            n_matter++;

            double x   = psi_field.get_x(i, j);
            double psi = psi_field.get_data(i, j);

            double vphys = std::fabs(get_v_raw(i, j, x, psi) * x * psi * psi);

            if (vphys > max_vphys) {
                max_vphys = vphys;
                bad_x = x;
                bad_y = psi_field.get_y(i, j);
            }
            if (vphys > 1.) n_bad++;
        }
    }

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  SANITY  |v_phys| max = " << max_vphys
              << std::setprecision(4)
              << "  (psi^2 v r, at x = " << bad_x << ", y = " << bad_y << ")\n";
    std::cout << std::defaultfloat;
    std::cout << "  SANITY  superluminal cells = " << n_bad
              << " of " << n_matter << " matter cells"
              << (n_bad ? "   *** NOT PHYSICAL ***" : "") << "\n\n";

    return n_bad == 0;
}


// -----------------------------------------------------------------------
// Save2 — same format as Grid::save_state()
// -

void Source::save(const std::string& out_dir, int ng,
                  const Field& psi_field) const {
    int nx = m_nxg - 2*ng;
    int ny = m_nyg - 2*ng;

    // the default stream buffer is a few KB and every << goes through
    // locale-aware formatting. a 4MB buffer and the C locale between them
    // roughly halve the time this takes, for identical output
    std::vector<std::vector<char>> bufs(6, std::vector<char>(1 << 22));
    std::ofstream f_rho_raw, f_rho_eff, f_src_psi, f_v_raw, f_v_eff, f_src_W;
    {
        std::ofstream* fs[6] = {&f_rho_raw, &f_rho_eff, &f_src_psi,
                                &f_v_raw,   &f_v_eff,   &f_src_W};
        const char* names[6] = {"/rho-raw.dat", "/rho-eff.dat", "/source-psi.dat",
                                "/v-raw.dat",   "/v-eff.dat",   "/source-W.dat"};
        for (int k = 0; k < 6; k++) {
            fs[k]->rdbuf()->pubsetbuf(bufs[k].data(), bufs[k].size());
            fs[k]->imbue(std::locale::classic());
            fs[k]->open(out_dir + names[k]);
        }
    }

    // set precision on all
    for (auto* f : {&f_rho_raw, &f_rho_eff, &f_src_psi,
                    &f_v_raw,   &f_v_eff,   &f_src_W})
        *f << std::scientific << std::setprecision(10);

    for (int j = ng; j < ng + ny; j++) {
        for (int i = ng; i < ng + nx; i++) {
            double x      = psi_field.get_x(i, j);
            double psi_ij = psi_field.get_data(i, j);

            // a char, not a std::string: this line ran once per cell, and
            // the heap allocation it used to do cost more than 3x the rest
            // of the write put together
            const char sep = (i < ng + nx - 1) ? ',' : '\n';

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

