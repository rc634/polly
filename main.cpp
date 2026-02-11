#include <iostream>
#include "field.hpp"
#include "grid.hpp"
#include "multigrid.hpp"
#include "problem.hpp"
#include <cmath>


int main() {

    Multigrid mg;

    double old_delta;
    double new_delta;
    double old_integral;

    mg.init();
    mg.hello();

    // v cycle controlls 
    int i_top = mg.i_fine;
    int i_bot = mg.i_coarse;
    int n_trigger = 200;

    // set zeros for the solver to start from
    mg.initial_data();
    mg.save_data("g");

    

    // for (int n=1; n<=20000; n++) {
    //     // lets goooooo!

    //     mg.v_cycle(i_top,i_bot);
    //     //mg.w_cycle();

    //     new_delta = mg.fine_grid_ptr->m_delta;

    //     if (n%100==0)
    //     {
    //         std::cout << ">- VCYCLE " << n << " : ";
    //         std::cout << "integral " << mg.fine_grid_ptr->field_integral();
    //         std::cout << " : delta " << mg.fine_grid_ptr->m_delta << "\n";
    //         std::cout << " : d-delta " << new_delta-old_delta << "\n";
    //     }

    //     if (n%n_trigger==0 && abs(new_delta-old_delta)<10e-18) {
    //         std::cout << "#################\n";
    //         std::cout << "# Dropping grid "<< i_bot << "!\n";
    //         std::cout << "#################\n";
    //         i_bot += 1;
    //         // stop v-cycling if we have one level left
    //         if (i_top == i_bot) break;
    //     }
        
    //     old_delta = new_delta;
    // }

    for (int n=1; n<=100; n++) {
        mg.w_cycle();
        std::cout << ">- WCYCLE " << n << " : ";
        std::cout << "integral W " << mg.fine_grid_ptr->int_W();
        std::cout << " : integral psi " << mg.fine_grid_ptr->int_psi() << "\n";
        // std::cout << " : delta " << mg.fine_grid_ptr->m_delta << "\n";
        std::cout << " : delta integral " 
                << mg.fine_grid_ptr->field_integral()-old_integral << "\n";

        if (abs(mg.fine_grid_ptr->field_integral()-old_integral)<10e-20) {
            break;
        }
        old_integral = mg.fine_grid_ptr->field_integral();
    }



    // for (int n=1; n<=100000; n++) {
    //     // lets goooooo!

    //     mg.refine();

    //     new_delta = mg.fine_grid_ptr->m_delta;

    //     if (n%10000==0)
    //     {
    //         std::cout << ">- REFINEMENT " << n << " : ";
    //         std::cout << "integral " << mg.fine_grid_ptr->field_integral();
    //         std::cout << " : delta " << mg.fine_grid_ptr->m_delta << "\n";
    //         std::cout << " : d-delta " << new_delta-old_delta << "\n";
    //     }

    //     old_delta = new_delta;
    // }



    for (int n=1; n<=600000; n++) {
        // lets goooooo!

        mg.refine();

        if (n%10000==0) {
            std::cout << ">- REFINEMENT " << n << " : ";
            std::cout << "integral W " << mg.fine_grid_ptr->int_W();
            std::cout << " : integral psi " << mg.fine_grid_ptr->int_psi() << "\n";
            // std::cout << " : delta " << mg.fine_grid_ptr->m_delta << "\n";
            std::cout << " : delta integral " 
                    << mg.fine_grid_ptr->field_integral()-old_integral << "\n";

            if (abs(mg.fine_grid_ptr->field_integral()-old_integral)<10e-20) {
                break;
        }
        }
        old_integral = mg.fine_grid_ptr->field_integral();
    }




    mg.save_data("f");

    mg.fine_grid_ptr->save_state();

    return 0;
}

