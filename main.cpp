#include <iostream>
#include "field.hpp"
#include "grid.hpp"
#include "multigrid.hpp"
#include "problem.hpp"

int main() {

    Multigrid mg;

    mg.init();

    mg.hello();

    mg.v_cycle();
    mg.v_cycle();

    mg.save_data();

    return 0;
}

