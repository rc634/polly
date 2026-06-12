#include <iostream>
#include <filesystem>
#include <cmath>
#include "multigrid.hpp"
#include "params.hpp"

struct RunConfig {
    double      rho_0;
    double      v_0;
    double      a;      // source height
    double      b;      // source width
    std::string name;
};

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
    mg.up_cycle(mg.p.tol);
    mg.fine_grid_ptr->save_state();
}

int main() {

    // -------------------------------------------------------
    // Parameter sweep — edit this list to change the batch
    // -------------------------------------------------------
    std::vector<RunConfig> runs = {
        //  rho_0   v_0   a     b          name
        {   3.02,   0.0,  1.0,  1.0,  "batch/rho3-02" },
        {   3.04,   0.0,  1.0,  1.0,  "batch/rho3-04" },
        {   3.06,   0.0,  1.0,  1.0,  "batch/rho3-06" },
        {   3.08,   0.0,  1.0,  1.0,  "batch/rho3-08" },
        {   3.10,   0.0,  1.0,  1.0,  "batch/rho3-10" },
        {   3.12,   0.0,  1.0,  1.0,  "batch/rho3-12" },
        {   3.14,   0.0,  1.0,  1.0,  "batch/rho3-14" },
        {   3.16,   0.0,  1.0,  1.0,  "batch/rho3-16" },
        {   3.18,   0.0,  1.0,  1.0,  "batch/rho3-18" },
    };
    // -------------------------------------------------------

    for (auto& cfg : runs)
        run(cfg);

    std::cout << "\nAll runs complete.\n";
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
