#include <iostream>
#include "field.hpp"
#include "grid.hpp"
#include "multigrid.hpp"
#include "problem.hpp"

int main() {

    Multigrid mg;

    double old_L2;
    double new_L2;

    mg.init();

    mg.hello();

    for (int n=1; n<=2000; n++) {
        // lets goooooo!

        mg.v_cycle();

        new_L2 = mg.fine_grid_ptr->L2norm();

        if (n%100==0)
        {
            std::cout << ">- VCYCLE " << n << " : ";
            //std::cout << "L2 norm " << new_L2;
            std::cout << "integral " << mg.fine_grid_ptr->field_integral();
            std::cout << " : delta " << mg.fine_grid_ptr->m_delta << "\n";
            // std::cout << " : remainder " << new_L2-old_L2 << "\n";
        }

        
        
        old_L2 = new_L2;
    }

    for (int n=1; n<=1000; n++) {
        // lets goooooo!

        mg.refine();

        new_L2 = mg.fine_grid_ptr->L2norm();

        if (n%100==0)
        {
            std::cout << ">- VCYCLE " << n << " : ";
            //std::cout << "L2 norm " << new_L2;
            std::cout << "integral " << mg.fine_grid_ptr->field_integral();
            std::cout << " : delta " << mg.fine_grid_ptr->m_delta << "\n";
            // std::cout << " : remainder " << new_L2-old_L2 << "\n";
        }

        old_L2 = new_L2;
    }

    // hacked override data with initial data for plotting 
    // mg.initial_data();

    mg.save_data();

    return 0;
}

