#include <fstream>
#include <string>
#include "grid.hpp"
#include <iomanip> // for precision output
#include <algorithm> // for std::min

// for pi 
#define _USE_MATH_DEFINES
#include <cmath>


Grid::Grid() {
}

void Grid::init(int i) {
    // numerical gridpoint stuff 
    grid_n = (int) i;
    nx = (int) pow(2,grid_n) * p.nx;
    ny = (int) pow(2,grid_n) * p.ny;
    // cell centered
    dx = (p.xU-p.xL)/nx;
    dy = (p.yU-p.yL)/ny;

    // deal with ghosts 
    ng = p.ng;
    nxg = nx + 2*ng;
    nyg = ny + 2*ng;
    // flattened array index
    n_flat = nxg*nyg;

    // indices 
    imin = ng;
    imax = ng + nx;
    jmin = ng;
    jmax = ng + ny;

    // feild data stuff 
    W.init(grid_n);
    psi.init(grid_n);

    // courant friendly timestep
    m_dt = p.CFL*std::min(dx,dy)*std::min(dx,dy);
    m_SOR = 1.;
    
    // verbose 
    std::cout << " - grid " << grid_n << " - dx = " << dx << "\n";
}

void Grid::ID_gaussian(double x0, double y0, double sig) {
    for (int j = jmin; j < jmax; j++) {
        for (int i = imin; i < imax; i++) {
    // for (int j = 0; j < nyg; j++) {
    //     for (int i = 0; i < nxg; i++) {
            // r^2 from centre of gaussian
            double rr = pow((x0-W.get_x(i,j)),2) + pow((y0-W.get_y(i,j)),2);
            // gaussian
            double val = exp(-(rr)/(2.*sig*sig));
            // set data 
            W.set_data(val,i,j);
        }
    }
}

void Grid::ID_pdisc(double R, double density) {
    // for (int j = jmin; j < jmax; j++) {
    //     for (int i = imin; i < imax; i++) {
    for (int j = 0; j < nyg; j++) {
        for (int i = 0; i < nxg; i++) {
            // r^2 from centre of gaussian
            double x = W.get_x(i,j);
            double y = W.get_y(i,j);
            double rr = pow((0.0-x),2) + pow((0.5-y),2);
            double out = 0.0;
            // disc of densiy for radii < R
            if (rr < R*R) {
                out = density*(rr/6. - 0.5*R*R);
                W.set_data(out,i,j);
            } 
            else {
                out = density*( - R*R*R / (3.*sqrt(rr)) );
                W.set_data(out,i,j);
            }
        }
    }
}

void Grid::ID_zeros() {

    //////////////////////////////////////////
    // // excludes ghosts 
    // for (int j = jmin; j < jmax; j++) {
    //     for (int i = imin; i < imax; i++) {

    //////////////////////////////////////////
    // includes ghosts 
    for (int j = 0; j < nyg; j++) {
        for (int i = 0; i < nxg; i++) {

        // zero initial data 
            W.set_data(0.,i,j);
            psi.set_data(1.,i,j);
        }
    }
}

double Grid::field_integral() {
    // can add metrics if wanted 
    double integral = 0.;
    double root_g = 1.;
    for (int j = jmin; j < jmax; j++) {
        for (int i = imin; i < imax; i++) {
            // root_g should = x
            integral += W.get_data(i,j) * root_g * dx * dy;
            integral += psi.get_data(i,j) * root_g * dx * dy;
        }
    }
    return integral;
}

double Grid::int_W() {
    // can add metrics if wanted 
    double integral = 0.;
    double root_g = 1.;
    for (int j = jmin; j < jmax; j++) {
        for (int i = imin; i < imax; i++) {
            // root_g should = x
            integral += W.get_data(i,j) * root_g * dx * dy;
        }
    }
    return integral;
}

double Grid::int_psi() {
    // can add metrics if wanted 
    double integral = 0.;
    double root_g = 1.;
    for (int j = jmin; j < jmax; j++) {
        for (int i = imin; i < imax; i++) {
            // root_g should = x
            integral += psi.get_data(i,j) * root_g * dx * dy;
        }
    }
    return integral;
}


void Grid::relax() {
    // One Relaxation step 
    
    // loop over all live cells 
    for (int j = jmin; j < jmax; j++) {
        for (int i = imin; i < imax; i++) {
            // time evolution step

            // const 
            double pi = M_PI;

            // coordinates 
            double x = W.get_x(i,j);
            double y = W.get_y(i,j);
            double r = sqrt(x*x + y*y);

            // field values 
            double Wij = W.get_data(i,j);
            double psij = psi.get_data(i,j);
            double dWdx = W.d1x(i,j);
            double dWdy = W.d1y(i,j);

            // source stats -- const elipsoid
            double height = 0.4;
            double width = 2.0; 
            // fractional coords wrt ellipsoid
            double xf = x/width;
            double yf = y/height;
            double rf = sqrt(xf*xf + yf*yf);

            // source terms
            double src_psi = 0.; // psi source
            double src_W = 0.; // W source
            double rho = 4.; // density source 
            double omega = 1.; // rotation source
            double packet = 0.; // overall shape of source 

            // if inside ellipse
            if (xf*xf + yf*yf < 1.) {
                packet = 1 + cos(rf * pi);
                packet *= packet;
                src_W = 8. * pi * x * rho * omega * packet;
                src_psi = - 2. * pi * rho * packet;
            }
            src_W += Wij/x/x;
            src_psi += - 0.5 * pow(psij,-7) * (
                dWdy*dWdy + dWdx*dWdx 
                - 2.*dWdx*Wij/x + Wij*Wij/x/x
            );

            // poisson eqaution
            double dW_dt = W.cylindrical_laplacian(i,j,x) - src_W;
            double dpsi_dt = psi.cylindrical_laplacian(i,j,x) - src_psi;

            // timestep f_new = f_old + dt * df/dt
            double new_W = W.get_data(i,j) + dW_dt * m_dt * m_SOR;
            double new_psi = psi.get_data(i,j) + dpsi_dt * m_dt * m_SOR;

            // set new-data stage 
            W.set_new_data(new_W,i,j);
            psi.set_new_data(new_psi,i,j);
        }
    }

    // debugging
    // W.copy_old_ghosts_to_new();


    // immediately after creating new data must fill its ghosts!
    // before diffing!
    fill_all_ghosts();

    // change in new solution, L2 norm
    m_delta = W.delta_data() + psi.delta_data();

    // Save over the old data - 
    // - with the finished new data
    W.save_new_data();
    psi.save_new_data();

    fill_all_ghosts();
}

// newer mixed type boundary conditions based on asymptotic expansion
void Grid::fill_all_ghosts() {
    // two ghosts thickness hard coded!

    // warning only one thickness of ghost computed for outer!
    // don't use 4th order stencils unless we fix that

    // loop left boundary -- symmetric psi, antisymmetric W
    for (int j = 0; j < nyg; j++) {
        W.set_data(-W.get_data(3,j),0,j);
        W.set_data(-W.get_data(2,j),1,j);
        W.set_new_data(-W.get_new_data(3,j),0,j);
        W.set_new_data(-W.get_new_data(2,j),1,j);
        psi.set_data(psi.get_data(3,j),0,j);
        psi.set_data(psi.get_data(2,j),1,j);
        psi.set_new_data(psi.get_new_data(3,j),0,j);
        psi.set_new_data(psi.get_new_data(2,j),1,j);
    }

    // loop bottom boundary -- symmetric psi and W
    for (int i = 0; i < nxg; i++) {
        W.set_data(W.get_data(i,2), i, 1);
        W.set_data(W.get_data(i,3), i, 0);
        W.set_new_data(W.get_new_data(i,2), i, 1);
        W.set_new_data(W.get_new_data(i,3), i, 0);
        psi.set_data(psi.get_data(i,2), i, 1);
        psi.set_data(psi.get_data(i,3), i, 0);
        psi.set_new_data(psi.get_new_data(i,2), i, 1);
        psi.set_new_data(psi.get_new_data(i,3), i, 0);
    }

    // outer boundaries

    // psi ~ 1 + M/2R   R is spherical radius

    // W ~ J sin theta / R^2     where J is some kind of rotation magnitude constant

    // declarations for mixed/robin bonudaries 
    double x2 = 0.; // cylinder radius 
    double y2 = 0.; // cylinder z
    double x3 = 0.; // cylinder radius 
    double y3 = 0.; // cylinder z
    // 3 is inner ghost, 4 is outer ghost
    double psi1 = 1., psi2 = 1., psi3 = 1., psi4 = 1.; 
    double W1=0., W2=0., W3=0., W4=0.; 
    int i4 = nxg-1, i3 = nxg-2, i2 = nxg-3, i1 = nxg-4;
    int j4 = nyg-1, j3 = nyg-2, j2 = nyg-3, j1 = nyg-4;

    // loop right boundary
    //  - x const, partial x calculated
    for (int j = 0; j < nyg; j++) {
        // calculate values for mixed deriv conditions
        x2 = W.get_x(i2,j);
        y2 = W.get_y(i2,j);
        x3 = W.get_x(i3,j);
        y3 = W.get_y(i3,j);

        psi1 = psi.get_data(i1,j); 
        psi2 = psi.get_data(i2,j);
        psi3 = psi1 - 2. * dx * x2 * (psi2-1.) / (x2*x2 + y2*y2);
        psi4 = psi2 - 2. * dx * x3 * (psi3-1.) / (x3*x3 + y3*y3);

        W1 = W.get_data(i1,j); 
        W2 = W.get_data(i2,j);
        W3 = W1 + 2. * dx * (y2*y2 - 2.*x2*x2) * W2 / (y2*y2*x2 + x2*x2*x2);
        W4 = W2 + 2. * dx * (y3*y3 - 2.*x3*x3) * W3 / (y3*y3*x3 + x3*x3*x3);

        // inner ghost
        W.set_data(W3, i3, j);
        W.set_new_data(W3, i3, j);
        psi.set_data(psi3, i3, j);
        psi.set_new_data(psi3, i3, j);
        // outer ghost 
        W.set_data(W4, i4, j);
        W.set_new_data(W4, i4, j);
        psi.set_data(psi4, i4, j);
        psi.set_new_data(psi4, i4, j);
    }

    // loop top boundary
    //  - y const, partial y calculated
    for (int i = 0; i < nxg; i++) {
        // calculate values for mixed deriv conditions
        x2 = W.get_x(i,j2);
        y2 = W.get_y(i,j2);
        x3 = W.get_x(i,j3);
        y3 = W.get_y(i,j3);

        psi1 = psi.get_data(i,j1); 
        psi2 = psi.get_data(i,j2);
        psi3 = psi1 - 2. * dx * x2 * (psi2-1.) / (x2*x2 + y2*y2);
        psi4 = psi2 - 2. * dx * x3 * (psi3-1.) / (x3*x3 + y3*y3);

        W1 = W.get_data(i,j1); 
        W2 = W.get_data(i,j2);
        W3 = W1 - 6. * dy * y2 * W2 / (y2*y2 + x2*x2);
        W4 = W2 - 6. * dy * y3 * W3 / (y3*y3 + x3*x3);

        // inner ghost 
        W.set_data(W3, i, j3);
        W.set_new_data(W3, i, j3);
        psi.set_data(psi3, i, j3);
        psi.set_new_data(psi3, i, j3);
        // outer ghost 
        W.set_data(W4, i, j4);
        W.set_new_data(W4, i, j4);
        psi.set_data(psi4, i, j4);
        psi.set_new_data(psi4, i, j4);
    }

}

void Grid::fill_all_ghosts_dirichlet() {


    // loop left boundary -- symmetric
    for (int j = 0; j < nyg; j++) {
        W.set_data(W.get_data(3,j),0,j);
        W.set_data(W.get_data(2,j),1,j);
        W.set_new_data(W.get_new_data(3,j),0,j);
        W.set_new_data(W.get_new_data(2,j),1,j);
        psi.set_data(psi.get_data(3,j),0,j);
        psi.set_data(psi.get_data(2,j),1,j);
        psi.set_new_data(psi.get_new_data(3,j),0,j);
        psi.set_new_data(psi.get_new_data(2,j),1,j);
    }

    // loop bottom boundary -- symetric 
    for (int i = 0; i < nxg; i++) {
        W.set_data(W.get_data(i,2), i, 1);
        W.set_data(W.get_data(i,3), i, 0);
        W.set_new_data(W.get_new_data(i,2), i, 1);
        W.set_new_data(W.get_new_data(i,3), i, 0);
        psi.set_data(psi.get_data(i,2), i, 1);
        psi.set_data(psi.get_data(i,3), i, 0);
        psi.set_new_data(psi.get_new_data(i,2), i, 1);
        psi.set_new_data(psi.get_new_data(i,3), i, 0);
    }

    // loop right boundary
    for (int j = 0; j < nyg; j++) {
        W.set_data(0., nxg-2, j);
        W.set_data(0., nxg-1, j);
        W.set_new_data(0., nxg-2, j);
        W.set_new_data(0., nxg-1, j);
        psi.set_data(1., nxg-2, j);
        psi.set_data(1., nxg-1, j);
        psi.set_new_data(1., nxg-2, j);
        psi.set_new_data(1., nxg-1, j);
    }

    // loop top boundary
    for (int i = 0; i < nxg; i++) {
        W.set_data(0., i, nyg-2);
        W.set_data(0., i, nyg-1);
        W.set_new_data(0., i, nyg-2);
        W.set_new_data(0., i, nyg-1);
        psi.set_data(1., i, nyg-2);
        psi.set_data(1., i, nyg-1);
        psi.set_new_data(1., i, nyg-2);
        psi.set_new_data(1., i, nyg-1);
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

    // set savae resolution
    file << std::scientific << std::setprecision(p.save_precision);

    // auto f_ = W; // lazy 

    // write 
    // excluding ghost cells 
    for (int j = jmin; j < jmax; j++) {
        for (int i = imin; i < imax; i++) {
    // // including ghost cells 
    // for (int j = 0; j < nyg; j++) {
    //     for (int i = 0; i < nxg; i++) {
            file << W.get_data(i,j) ;
            //if (i==nxg-1) {
            if (i==imax-1) {
                break;
            }
            file << ",";
        }
        file << "\n";
    }

    file.close();
}

void Grid::save_state() {
    std::string path1 = "data/W.dat";
    std::ofstream file1(path1);
    std::string path2 = "data/psi.dat";
    std::ofstream file2(path2);

    if (!file1) {
        std::cerr << "Error: could not open file " << path1 << "\n";
        return;
    }

    if (!file2) {
        std::cerr << "Error: could not open file " << path2 << "\n";
        return;
    }

    std::cout << "Saving " << path1 << std::endl;
    std::cout << "Saving " << path2 << std::endl;

    // set save resolution
    file1 << std::scientific << std::setprecision(p.save_precision);
    file2 << std::scientific << std::setprecision(p.save_precision);

    // write 
    // excluding ghost cells 
    for (int j = jmin; j < jmax; j++) {
        for (int i = imin; i < imax; i++) {
    // // including ghost cells 
    // for (int j = 0; j < nyg; j++) {
    //     for (int i = 0; i < nxg; i++) {
            file1 << W.get_data(i,j) ;
            file2 << psi.get_data(i,j) ;
            //if (i==nxg-1) {
            if (i==imax-1) {
                break;
            }
            file1 << ",";
            file2 << ",";
        }
        file1 << "\n";
        file2 << "\n";
    }

    file1.close();
    file2.close();
}




void Grid::hello() {
    std::cout << "Hello from Grid :  " << grid_n << "\n"; 
    W.hello();
}
