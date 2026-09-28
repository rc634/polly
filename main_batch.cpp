#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <filesystem>
#include <cmath>
#include <chrono>
#include "multigrid.hpp"
#include "params.hpp"

struct RunConfig {
    double      rho_0;
    double      v_0;
    double      a;      // source height
    double      b;      // source width
    std::string name;
};

// -----------------------------------------------------------------------
// Dumps the run parameters next to the grid data so downstream tools
// (cauchy-tools) know the source shape without re-deriving it.
// Format is one "key value" pair per line, so a reader is just
//   while (f >> key >> value) ...
// -----------------------------------------------------------------------
static void save_run_params(const RunConfig& cfg, const Grid& g) {
    std::string path = cfg.name + "/params.dat";
    std::ofstream f(path);

    if (!f) {
        std::cerr << "Error: could not open file " << path << "\n";
        return;
    }

    std::cout << "Saving " << path << std::endl;

    // aspect ratio of the short axis to the long one, and the
    // eccentricity of the spheroid it describes
    double q   = std::min(cfg.a, cfg.b) / std::max(cfg.a, cfg.b);
    double ecc = sqrt(1. - q*q);

    f << std::scientific << std::setprecision(g.p.save_precision);

    // source
    f << "rho_0 " << cfg.rho_0 << "\n";
    f << "v_0 "   << cfg.v_0   << "\n";
    f << "a "     << cfg.a     << "\n";   // height, y semi-axis
    f << "b "     << cfg.b     << "\n";   // width,  x semi-axis
    f << "ecc "   << ecc       << "\n";
    f << "M_ADM " << g.m_source.m_M_ADM << "\n";

    // grid the .dat files live on (finest level, ghosts excluded)
    f << "nx " << g.nx   << "\n";
    f << "ny " << g.ny   << "\n";
    f << "xL " << g.p.xL << "\n";
    f << "xU " << g.p.xU << "\n";
    f << "yL " << g.p.yL << "\n";
    f << "yU " << g.p.yU << "\n";
}

static void run(const RunConfig& cfg) {
    std::cout << "\n";
    std::cout << "*=========================================*\n";
    std::cout << "| Run : " << cfg.name << "\n";
    std::cout << "| rho_0 = " << cfg.rho_0 << "  v_0 = " << cfg.v_0 << "  a = " << cfg.a << "  b = " << cfg.b << "\n";
    std::cout << "*=========================================*\n";

    std::filesystem::create_directories(cfg.name);

    Multigrid mg;
    for (int i = 0; i < mg.p.num_grids; i++) {
        mg.grids[i].p.out_dir = cfg.name;
        mg.grids[i].m_source.set_source_params(cfg.rho_0, cfg.v_0, cfg.a, cfg.b);
    }

    mg.init();
    mg.initial_data();
    mg.up_cycle(mg.p.res_tol);
    mg.fine_grid_ptr->save_state();
    save_run_params(cfg, *mg.fine_grid_ptr);
}

int main() {

    auto t_start = std::chrono::steady_clock::now();

    // -------------------------------------------------------
    // Parameter sweep — edit this list to change the batch
    // -------------------------------------------------------
    const int n = 1;
    double a_bag[n];
    double b_bag[n];
    double a_pan[n];
    double b_pan[n];
    double ecc_i[n] = {0.99};
    double spins[7] = {0.1, 0.3, 0.6, 1.0, 2.0, 4.0, 8.0};
    for (size_t i = 0; i < n; i++)
    {
        // height a, widths b, b
        double a_over_b = sqrt(1.-ecc_i[i]*ecc_i[i]);

        // constant volume ~ a*b*b=1

        // pancake has a < b
        b_pan[i] = pow(a_over_b,-1./3.);
        a_pan[i] = pow(a_over_b, 2./3.);

        // baguette has a > b
        b_bag[i] = pow(a_over_b, 1./3.);
        a_bag[i] = pow(a_over_b,-2./3.);
    }

    // e99 bag about 7.6
    // e99 pan about 5.2
    // e99 sph about 3.2
    // const double rho_bag = 5.;
    // const double rho_pan = 4.;
    // const double rho_bag99 = 7.7;
    // const double rho_pan99 = 5.3;
    // const double rho_sph = 3.1;

    std::vector<RunConfig> runs = {
        // rho_0  v_0   a     b          name
        // {   3.08,   0.0,  1.0,  1.0,  "batch/paper-sph/rho3-08" }, 
        // {   3.10,   0.0,  1.0,  1.0,  "batch/paper-sph/rho3-10" },
        // {   3.12,   0.0,  1.0,  1.0,  "batch/paper-sph/rho3-12" }, 
        // {   3.14,   0.0,  1.0,  1.0,  "batch/paper-sph/rho3-14" }, 
        // {   3.16,   0.0,  1.0,  1.0,  "batch/paper-sph/rho3-16" }, 
        // {   3.18,   0.0,  1.0,  1.0,  "batch/paper-sph/rho3-18" }, 
        // {   3.20,   0.0,  1.0,  1.0,  "batch/paper-sph/rho3-20" }, //comfortable horizon
        {   5.0,   0.0,  1.0,  1.0,  "batch/paper-sph/rho5-0" }, 

        // {   7.0,   0.0,  a_bag[0],  b_bag[0],  "batch/paper-bag/rho7-00" }, // no
        // {   7.1,   0.0,  a_bag[0],  b_bag[0],  "batch/paper-bag/rho7-10" }, 
        // {   7.2,   0.0,  a_bag[0],  b_bag[0],  "batch/paper-bag/rho7-20" }, 
        // {   7.3,   0.0,  a_bag[0],  b_bag[0],  "batch/paper-bag/rho7-30" }, 
        // {   7.4,   0.0,  a_bag[0],  b_bag[0],  "batch/paper-bag/rho7-40" }, 
        // {   7.5,   0.0,  a_bag[0],  b_bag[0],  "batch/paper-bag/rho7-50" }, // yes


        // {   4.5,   0.0,  a_pan[0],  b_pan[0],  "batch/paper-pan/rho4-50" }, //no
        // {   4.6,   0.0,  a_pan[0],  b_pan[0],  "batch/paper-pan/rho4-60" }, 
        // {   4.7,   0.0,  a_pan[0],  b_pan[0],  "batch/paper-pan/rho4-70" }, 
        // {   4.8,   0.0,  a_pan[0],  b_pan[0],  "batch/paper-pan/rho4-80" }, 
        // {   4.9,   0.0,  a_pan[0],  b_pan[0],  "batch/paper-pan/rho4-90" }, 
        // {   5.0,   0.0,  a_pan[0],  b_pan[0],  "batch/paper-pan/rho5-00" }, //
        
    };

    for (auto& cfg : runs)
        run(cfg);

    auto t_end = std::chrono::steady_clock::now();
    long long total_s = std::chrono::duration_cast<std::chrono::seconds>(t_end - t_start).count();
    long long mins = total_s / 60;
    long long secs = total_s % 60;

    std::cout << "\nAll runs complete.\n";
    std::cout << "Total runtime: " << mins << " min " << secs << " sec\n";
    return 0;
}




// std::vector<RunConfig> runs = {
//         //  rho_0  v_0   a     b          name
//         {   4.0,   0.0,  1.0,  1.0,  "batch/rho4-0" },
//         {   3.8,   0.0,  1.0,  1.0,  "batch/rho3-8" },
//         {   3.6,   0.0,  1.0,  1.0,  "batch/rho3-6" },
//         {   3.4,   0.0,  1.0,  1.0,  "batch/rho3-4" },
//         {   3.2,   0.0,  1.0,  1.0,  "batch/rho3-2" },
//         {   3.0,   0.0,  1.0,  1.0,  "batch/rho3-0" },
//         {   2.8,   0.0,  1.0,  1.0,  "batch/rho2-8" },
//         {   2.6,   0.0,  1.0,  1.0,  "batch/rho2-6" },
//         {   2.4,   0.0,  1.0,  1.0,  "batch/rho2-4" },
//         {   2.2,   0.0,  1.0,  1.0,  "batch/rho2-2" },
//         {   2.0,   0.0,  1.0,  1.0,  "batch/rho2-0" },
//     };






        // {   6.1,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho006-1" },
        // {   6.2,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho006-2" },
        // {   6.3,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho006-3" },
        // {   6.4,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho006-4" },
        // {   6.5,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho006-5" },
        // {   6.6,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho006-6" },
        // {   6.7,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho006-7" },
        // {   6.8,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho006-8" },
        // {   6.9,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho006-9" },
        // {   7.0,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho007-0" },
        // {   7.1,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho007-1" },
        // {   7.2,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho007-2" },
        // {   7.3,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho007-3" },
        // {   7.4,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho007-4" },
        // {   7.5,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho007-5" },
        // {   7.6,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho007-6" },
        // {   7.7,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho007-7" },
        // {   7.8,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho007-8" },
        // {   7.9,   0.0,  a_bag,  b_bag,  "batch/bagfine/rho007-9" },
        // // pancake 
        // {   5.3,   0.0,  a_pan,  b_pan,  "batch/panfine/rho005-3" },
        // {   5.4,   0.0,  a_pan,  b_pan,  "batch/panfine/rho005-4" },
        // {   5.5,   0.0,  a_pan,  b_pan,  "batch/panfine/rho005-5" },



// std::vector<RunConfig> runs = {
//         //  rho_0   v_0   a     b          name
//         // coarse density 
//         {   4.0,   0.0,  1.0,  1.0,  "batch/coarse/rho4-0" },
//         {   3.8,   0.0,  1.0,  1.0,  "batch/coarse/rho3-8" },
//         {   3.6,   0.0,  1.0,  1.0,  "batch/coarse/rho3-6" },
//         {   3.4,   0.0,  1.0,  1.0,  "batch/coarse/rho3-4" },
//         {   3.2,   0.0,  1.0,  1.0,  "batch/coarse/rho3-2" },
//         {   3.0,   0.0,  1.0,  1.0,  "batch/coarse/rho3-0" },
//         {   2.8,   0.0,  1.0,  1.0,  "batch/coarse/rho2-8" },
//         {   2.6,   0.0,  1.0,  1.0,  "batch/coarse/rho2-6" },
//         {   2.4,   0.0,  1.0,  1.0,  "batch/coarse/rho2-4" },
//         {   2.2,   0.0,  1.0,  1.0,  "batch/coarse/rho2-2" },
//         {   2.0,   0.0,  1.0,  1.0,  "batch/coarse/rho2-0" },
//         // fine density 
//         {   3.02,   0.0,  1.0,  1.0,  "batch/fine/rho3-02" },
//         {   3.04,   0.0,  1.0,  1.0,  "batch/fine/rho3-04" },
//         {   3.06,   0.0,  1.0,  1.0,  "batch/fine/rho3-06" },
//         {   3.08,   0.0,  1.0,  1.0,  "batch/fine/rho3-08" },
//         {   3.10,   0.0,  1.0,  1.0,  "batch/fine/rho3-10" },
//         {   3.12,   0.0,  1.0,  1.0,  "batch/fine/rho3-12" },
//         {   3.14,   0.0,  1.0,  1.0,  "batch/fine/rho3-14" },
//         {   3.16,   0.0,  1.0,  1.0,  "batch/fine/rho3-16" },
//         {   3.18,   0.0,  1.0,  1.0,  "batch/fine/rho3-18" },
//         // sub spin
//         {   3.04,   1.00,  1.0,  1.0,  "batch/subspin/v1" },
//         {   3.04,   2.00,  1.0,  1.0,  "batch/subspin/v2" },
//         {   3.04,   4.00,  1.0,  1.0,  "batch/subspin/v4" },
//         {   3.04,   8.00,  1.0,  1.0,  "batch/subspin/v8" },
//         {   3.04,   16.0,  1.0,  1.0,  "batch/subspin/v16" },
//         {   3.04,   32.0,  1.0,  1.0,  "batch/subspin/v32" },
//         {   3.04,   64.0,  1.0,  1.0,  "batch/subspin/v64" },
//         // super spin
//         {   3.08,   1.00,  1.0,  1.0,  "batch/supspin/v1" },
//         {   3.08,   2.00,  1.0,  1.0,  "batch/supspin/v2" },
//         {   3.08,   4.00,  1.0,  1.0,  "batch/supspin/v4" },
//         {   3.08,   8.00,  1.0,  1.0,  "batch/supspin/v8" },
//         {   3.08,   16.0,  1.0,  1.0,  "batch/supspin/v16" },
//         {   3.08,   32.0,  1.0,  1.0,  "batch/supspin/v32" },
//         {   3.08,   64.0,  1.0,  1.0,  "batch/supspin/v64" },
//     };


// std::vector<RunConfig> runs = {
//         //  rho_0   v_0   a     b          name
//         // coarse density 
//         {   3.0,   0.0,  a_bag,  b_bag,  "batch/bag/rho003" },
//         {   10.0,   0.0,  a_bag,  b_bag,  "batch/bag/rho010" },
//         {   30.0,   0.0,  a_bag,  b_bag,  "batch/bag/rho030" },
//         {   100.0,   0.0,  a_bag,  b_bag,  "batch/bag/rho100" },
//         {   2.0,   0.0,  a_pan,  b_pan,  "batch/pan/rho002" },
//         {   4.0,   0.0,  a_pan,  b_pan,  "batch/pan/rho004" },
//         {   6.0,   0.0,  a_pan,  b_pan,  "batch/pan/rho006" },
//         {   8.0,   0.0,  a_pan,  b_pan,  "batch/pan/rho008" },
//         {   12.0,   0.0,  a_pan,  b_pan,  "batch/pan/rho012" },
//         {   20.0,   0.0,  a_pan,  b_pan,  "batch/pan/rho020" },
//         {   36.0,   0.0,  a_pan,  b_pan,  "batch/pan/rho036" },
//     };




// std::vector<RunConfig> runs = {
//         //  rho_0   v_0   a     b          name
//         // coarse density 
//         // baguette
//         {   rho_bag,   0.0,  a_bag[0],  b_bag[0],  "batch/bageccfine/ecc0940" },
//         {   rho_bag,   0.0,  a_bag[1],  b_bag[1],  "batch/bageccfine/ecc0950" },
//         {   rho_bag,   0.0,  a_bag[2],  b_bag[2],  "batch/bageccfine/ecc0960" },
//         {   rho_bag,   0.0,  a_bag[3],  b_bag[3],  "batch/bageccfine/ecc0970" },
//         {   rho_bag,   0.0,  a_bag[4],  b_bag[4],  "batch/bageccfine/ecc0980" },
//         // spinnin e99 baguette
//         {   rho_bag99, spins[0],  a_bag[5],  b_bag[5],  "batch/bagecc99spin/v01" },
//         {   rho_bag99, spins[1],  a_bag[5],  b_bag[5],  "batch/bagecc99spin/v03" },
//         {   rho_bag99, spins[2],  a_bag[5],  b_bag[5],  "batch/bagecc99spin/v06" },
//         {   rho_bag99, spins[3],  a_bag[5],  b_bag[5],  "batch/bagecc99spin/v10" },
//         {   rho_bag99, spins[4],  a_bag[5],  b_bag[5],  "batch/bagecc99spin/v20" },
//         {   rho_bag99, spins[5],  a_bag[5],  b_bag[5],  "batch/bagecc99spin/v40" },
//         {   rho_bag99, spins[6],  a_bag[5],  b_bag[5],  "batch/bagecc99spin/v80" },
//         // pancake 
//         {   rho_pan,   0.0,  a_pan[0],  b_pan[0],  "batch/paneccfine/ecc0940" },
//         {   rho_pan,   0.0,  a_pan[1],  b_pan[1],  "batch/paneccfine/ecc0950" },
//         {   rho_pan,   0.0,  a_pan[2],  b_pan[2],  "batch/paneccfine/ecc0960" },
//         {   rho_pan,   0.0,  a_pan[3],  b_pan[3],  "batch/paneccfine/ecc0970" },
//         {   rho_pan,   0.0,  a_pan[4],  b_pan[4],  "batch/paneccfine/ecc0980" },
//         // spinnin e99 pancake
//         {   rho_pan99, spins[0],  a_pan[5],  b_pan[5],  "batch/panecc99spin/v01" },
//         {   rho_pan99, spins[1],  a_pan[5],  b_pan[5],  "batch/panecc99spin/v03" },
//         {   rho_pan99, spins[2],  a_pan[5],  b_pan[5],  "batch/panecc99spin/v06" },
//         {   rho_pan99, spins[3],  a_pan[5],  b_pan[5],  "batch/panecc99spin/v10" },
//         {   rho_pan99, spins[4],  a_pan[5],  b_pan[5],  "batch/panecc99spin/v20" },
//         {   rho_pan99, spins[5],  a_pan[5],  b_pan[5],  "batch/panecc99spin/v40" },
//         {   rho_pan99, spins[6],  a_pan[5],  b_pan[5],  "batch/panecc99spin/v80" },
//         // spinnin sphere
//         {   rho_sph, spins[0],  a_pan[6],  b_pan[6],  "batch/sphspin/v01" },
//         {   rho_sph, spins[1],  a_pan[6],  b_pan[6],  "batch/sphspin/v03" },
//         {   rho_sph, spins[2],  a_pan[6],  b_pan[6],  "batch/sphspin/v06" },
//         {   rho_sph, spins[3],  a_pan[6],  b_pan[6],  "batch/sphspin/v10" },
//         {   rho_sph, spins[4],  a_pan[6],  b_pan[6],  "batch/sphspin/v20" },
//         {   rho_sph, spins[5],  a_pan[6],  b_pan[6],  "batch/sphspin/v40" },
//         {   rho_sph, spins[6],  a_pan[6],  b_pan[6],  "batch/sphspin/v80" },

//     };