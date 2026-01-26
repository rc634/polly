#pragma once
#include <iostream>
#include "grid.hpp"

class Multigrid {
public:

    Grid g;  

    Multigrid() { 
        std::cout << "Hello from Multigrid!\n"; 
    }

    void vcycle() { 
        std::cout << "Performing a dummy V-cycle\n"; 
    }
};

