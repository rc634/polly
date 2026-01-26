#pragma once
#include <iostream>
#include "field.hpp"

class Grid {
public:

    Field f; 
    
    Grid() { 
        std::cout << "Hello from Grid!\n"; 
    }
};

