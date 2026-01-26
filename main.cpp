#include <iostream>
#include "field.hpp"
#include "grid.hpp"
#include "multigrid.hpp"
#include "problem.hpp"

int main() {

    Multigrid mg;

    mg.vcycle();

    std::cout << "Hello from main!\n";
    
    return 0;
}

