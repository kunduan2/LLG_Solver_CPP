#include "utils.h"
#include "ExternalField.h"
#include "ExchangeField.h"
#include "DelmLlg.h"
#include "LlgSolver.h"

#include <iostream>
#include <random>
#include <fstream>
using namespace std;



int main() {

    /* Define all the parameters and variables used to solve the LLG */
    
    // Exchange
    const double Aexch = .250; 

    // Parameters to define the time dependent  External field H(t)
    const double t_switch = 50.0;  

    // Physical parameters for LLG
    const double gamma_gyro = 1.0;
    const double alpha = 0.1;

    // Temperature / thermal noise
    const double kB = 1.0;      // Boltzmann constant (choose units)
    const double T = 00.00;      // Temperature
    const double Ms = 1.0;     // Saturation magnetization
    const double V = 1.0;      // Volume
    const double D = alpha*kB*T/(gamma_gyro*Ms*V); // Noise amplitude 

    // Loop details
    int Nx = 5, Ny= 5;
    int Nt = 50;
    const double t = 0.0;
    const double dt = 0.01;

    // Bext-loop
    double Bmax = 5.5;
    double dB = .20;

    // make a grid for a distribution of magnetic moments in 2D
    double* Mx = new double[Nx*Ny];
    double* My = new double[Nx*Ny];
    double* Mz = new double[Nx*Ny];

    // initiatlize magnetization
    double mx = 1.0;
    double my = 0.0;
    double mz = 1.0;
    
    normalize(mx, my, mz);

    for (int i=0; i<Ny; i++){
        for (int j=0; j<Nx; j++){
            int ind_ij = convert_idx_2dto1D(i, j, Nx, Ny);            
            Mx[ind_ij]= mx;
            My[ind_ij]= my;
            Mz[ind_ij]= mz;
            // cout << Mz[ind_ij] << " ";
        }
    }

   print2Darr(Mx, Nx, Ny);


    LlgSolver llg_solver;
    llg_solver.solve(
        // Physics
        50.0,   // t_switch
        1.0,    // Bext
        1.0,    // Aexch
        1.0,    // D
        0.1,   // alpha
        1.0,    // gamma_gyro

        // grid/time
        Nt,      // Nt
        Nx,     // Nx
        Ny,     // Ny
        0.0,    // t
        0.01,   // dt

        // arrays
        Mx, //nullptr, // Mx
        My, //nullptr, // My
        Mz //nullptr  // Mz
    );


    // delete heap memory allocation
    delete[] Mx;
    delete[] My;
    delete[] Mz;

    return 0;
}
