#include <fstream>
#include <string>
#include "grid.hpp"
#include <iomanip> // for precision output
#include <algorithm> // for std::min


Grid::Grid() {
}

void Grid::init(int i) {
    grid_n = i;
    f1.init(grid_n);
    // courant friendly timestep
    m_dt = p.CFL*std::min(f1.dx,f1.dy)*std::min(f1.dx,f1.dy);
}

void Grid::ID_gaussian(double x0, double y0, double sig) {
    for (int j = f1.jmin; j < f1.jmax; j++) {
        for (int i = f1.imin; i < f1.imax; i++) {
    // for (int j = 0; j < f1.nyg; j++) {
    //     for (int i = 0; i < f1.nxg; i++) {
            // r^2 from centre of gaussian
            double rr = pow((x0-f1.get_x(i,j)),2) + pow((y0-f1.get_y(i,j)),2);
            // gaussian
            double val = exp(-(rr)/(2.*sig*sig));
            // set data 
            f1.set_data(val,i,j);
        }
    }
}

void Grid::ID_pdisc(double R, double density) {
    // for (int j = f1.jmin; j < f1.jmax; j++) {
    //     for (int i = f1.imin; i < f1.imax; i++) {
    for (int j = 0; j < f1.nyg; j++) {
        for (int i = 0; i < f1.nxg; i++) {
            // r^2 from centre of gaussian
            double x = f1.get_x(i,j);
            double y = f1.get_y(i,j);
            double rr = pow((0.0-x),2) + pow((0.5-y),2);
            double out = 0.0;
            // disc of densiy for radii < R
            if (rr < R*R) {
                out = density*(rr/6. - 0.5*R*R);
                f1.set_data(out,i,j);
            } 
            else {
                out = density*( - R*R*R / (3.*sqrt(rr)) );
                f1.set_data(out,i,j);
            }
        }
    }
}

void Grid::ID_zeros() {
    for (int j = f1.jmin; j < f1.jmax; j++) {
        for (int i = f1.imin; i < f1.imax; i++) {
    // for (int j = 0; j < f1.nyg; j++) {
    //     for (int i = 0; i < f1.nxg; i++) {
        // zero initial data 
            f1.set_data(0,i,j);
        }
    }
}

double Grid::L2norm() {
    double L2 = 0.;
    for (int k = 0; k < f1.n_flat; k++) {
        L2 += f1.data[k]*f1.data[k];
    }
    return sqrt(L2);
}

double Grid::field_integral() {
    // can add metrics if wanted 
    double integral = 0.;
    double root_g = 1.;
    for (int j = f1.jmin; j < f1.jmax; j++) {
        for (int i = f1.imin; i < f1.imax; i++) {
            // root g ...
            integral += f1.get_data(i,j) * root_g * f1.dx * f1.dy;
        }
    }
    return integral;
}


void Grid::relax() {
    // One Relaxation step 
    
    // loop over all live cells 
    for (int j = f1.jmin; j < f1.jmax; j++) {
        for (int i = f1.imin; i < f1.imax; i++) {
            // time evolution step

            // coordinate 
            double x = f1.get_x(i,j);
            double y = f1.get_y(i,j);
            double r = sqrt(x*x + y*y);

            // gaussian stats
            double x0 = 0.0;
            double y0 = 0.5;
            double r0 = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));
            double rmax = 0.1;
            double sigx = 0.1;
            double sigy = 0.15;

            // source terms
            double src_f1 = 0;
            if (r0 < rmax) {
                // src_f1 = exp(-(x-x0)*(x-x0)/(2*sigx*sigx))
                //        * exp(-(y-y0)*(y-y0)/(2*sigy*sigy));
                src_f1 = 100.;
            }

            // poisson eqaution
            // double df1_dt = f1.cartesian_laplacian(i,j) - src_f1;
            double df1_dt = f1.cylindrical_laplacian(i,j,x) - src_f1;

            // timestep f_new = f_old + dt * df/dt
            double new_f1 = f1.get_data(i,j) + df1_dt * m_dt;

            // set new-data stage 
            f1.set_new_data(new_f1,i,j);
        }
    }

    // debugging
    f1.copy_old_ghosts_to_new();


    // immediately after creating new data must fill its ghosts!
    // before diffing!
    // fill_all_ghosts();

    // change in new solution, L2 norm
    m_delta = f1.delta_data();

    // Save over the old data - 
    // - with the finished new data
    f1.save_new_data();

    // fill_all_ghosts();
}


void Grid::fill_all_ghosts() {
    // two ghosts thickness hard coded!

    // loop left boundary
    for (int j = 0; j < f1.nyg; j++) {
        f1.set_data(f1.get_data(3,j),0,j);
        f1.set_data(f1.get_data(2,j),1,j);
        f1.set_new_data(f1.get_new_data(3,j),0,j);
        f1.set_new_data(f1.get_new_data(2,j),1,j);
    }
    // fake 1/r boundary stuff 
    double dx = f1.dx, dy = f1.dy, dr=0.;
    double r = 1., x=0., y=0.;
    double nxg=f1.nxg, nyg=f1.nyg;
    // loop right boundary
    for (int j = 0; j < f1.nyg; j++) {
        y = -0.5+ j/f1.nyg;
        r = sqrt(1.+ y*y);
        dr = dx * abs(x) / r;
        f1.set_data( f1.get_data(nxg-3,j) * (1-dr/r), nxg-2, j);
        f1.set_data( f1.get_data(nxg-2,j) * (1-dr/(r+dr)), nxg-1, j);
        f1.set_new_data( f1.get_new_data(nxg-3,j) * (1-dr/r), nxg-2, j);
        f1.set_new_data( f1.get_new_data(nxg-2,j) * (1-dr/(r+dr)), nxg-1, j);
    }
    // loop top boundary
    for (int i = 0; i < f1.nxg; i++) {
        x = i/f1.nxg;
        y = 0.5;
        r = sqrt(0.25 + x*x);
        dr = dy * abs(y) / r;
        f1.set_data( f1.get_data(i,nyg-3) * (1-dr/r), i, nyg-2);
        f1.set_data( f1.get_data(i,nyg-2) * (1-dr/(r+dr)), i, nyg-1);
        f1.set_new_data( f1.get_new_data(i,nyg-3) * (1-dr/r), i, nyg-2);
        f1.set_new_data( f1.get_new_data(i,nyg-2) * (1-dr/(r+dr)), i, nyg-1);
    }
    // loop bottom boundary
    for (int i = 0; i < f1.nxg; i++) {
        x = i/f1.nxg;
        y = -0.5;
        r = sqrt(0.25 + x*x);
        dr = dy * abs(y) / r;
        f1.set_data( f1.get_data(i,2) * (1-dr/r), i,1);
        f1.set_data( f1.get_data(i,1) * (1-dr/(r+dr)), i,0);
        f1.set_new_data( f1.get_new_data(i,2) * (1-dr/r), i,1);
        f1.set_new_data( f1.get_new_data(i,1) * (1-dr/(r+dr)), i,0);
    }
}


void Grid::save_data(const std::string& filename) {
    std::string path = "data/" + filename + ".dat";
    std::ofstream file(path);

    if (!file) {
        std::cerr << "Error: could not open file " << path << "\n";
        return;
    }

    std::cout << "Saving " << path << std::endl;

    // 8 significant figures (not 8 decimal places!)
    file << std::scientific << std::setprecision(p.save_precision);

    // auto f_ = f1; // lazy 

    // write everything 
    // including ghost cells for now 
    for (int j = 0; j < f1.nyg; j++) {
        for (int i = 0; i < f1.nxg; i++) {
            file << f1.get_data(i,j) ;
            if (i==f1.nxg-1) {
                break;
            }
            file << ",";
        }
        file << "\n";
    }

    file.close();
}

void Grid::hello() {
    std::cout << "Hello from Grid :  " << grid_n << "\n"; 
    f1.hello();
}
