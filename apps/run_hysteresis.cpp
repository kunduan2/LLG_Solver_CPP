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
    const double Aexch = 0.0; 

    // Parameters to define the time dependent  External field H(t)
    const double t_switch = 50.0;  

    // Physical parameters for LLG
    const double gamma_gyro = .10;
    const double alpha = 0.2;

    // Temperature / thermal noise
    const double kB = 1.0;      // Boltzmann constant (choose units)
    const double T = 300.00;      // Temperature
    const double Ms = 1.0;     // Saturation magnetization
    const double V = 1.0;      // Volume
    const double D = alpha*kB*T/(gamma_gyro*Ms*V); // Noise amplitude 

    // Loop details
    int Nx = 50, Ny= 50;
    int Nt = 50;
    const double t = 0.0;
    const double dt = 0.01;

    // Bext-loop
    double Bmax = 10.2;
    double dB = .50;

    // Make a grid for a distribution of magnetic moments in 2D
    double* Mx = new double[Nx*Ny];
    double* My = new double[Nx*Ny];
    double* Mz = new double[Nx*Ny];

    // To track the loop
    int max_counter = 5*(Bmax/dB);


    // Open a file
    ofstream fout("test.dat");
    if (!fout.is_open()) {
        cerr << "Error opening file test.dat!" << endl;
        return 1; // or handle error appropriately
    }

    // initiatlize magnetization
    for (int i=0; i<Ny; i++){
        for (int j=0; j<Nx; j++){
            int ind_ij = convert_idx_2dto1D(i, j, Nx, Ny);            
            Mx[ind_ij]= 0.0;
            My[ind_ij]= 0.001;
            Mz[ind_ij]= 1;
            // cout << Mz[ind_ij] << " ";
        }
    }

    // save mean vals of initiatlize magnetization components
    fout << 0 << " " << mean(Mx, Nx*Ny) << " " << mean(My, Nx*Ny)<< " " << mean(Mz, Nx*Ny)<<"\n";   

    // define object
    LlgSolver llg_solver;

    // define variavles to be used in all for loops
    double mx = 1.0;
    double my = 0.0;
    double mz = 1.0;
    double Bext = 0.0;

    // 1st part of the Hysteresis loop
    for (int l=0; l<= (Bmax/dB); l++){
        Bext = l*dB;          
        llg_solver.solve(
            // Physics
            50.0,   // t_switch
            Bext,    // Bext
            Aexch,    // Aexch
            D,    // D
            alpha,   // alpha
            gamma_gyro,    // gamma_gyro

            // grid/time
            Nt,      // Nt
            Nx,     // Nx
            Ny,     // Ny
            t,    // t
            dt,   // dt

            // arrays
            Mx, //nullptr, // Mx
            My, //nullptr, // My
            Mz //nullptr  // Mz
        );
        mx = mean(Mx, Nx*Ny);
        my = mean(My, Nx*Ny);
        mz = mean(Mz, Nx*Ny);

        normalize(mx, my, mz);

        fout << Bext << " " << mx << " " << my << " " << mz<<"\n"; 
        // cout << Bext << " " << mx << " " << my << " " << mz<<"\n";  
        cout<< max_counter-- <<"\n";
    }

     // 2nd part of the Hysteresis loop
    for (int l=(Bmax/dB); l>= -(Bmax/dB); l--){  // for (double Bext = Bmax; Bext >= -Bmax; Bext -= dB){
        Bext = l*dB;  
                llg_solver.solve(
            // Physics
            50.0,   // t_switch
            Bext,    // Bext
            Aexch,    // Aexch
            D,    // D
            alpha,   // alpha
            gamma_gyro,    // gamma_gyro

            // grid/time
            Nt,      // Nt
            Nx,     // Nx
            Ny,     // Ny
            t,    // t
            dt,   // dt

            // arrays
            Mx, //nullptr, // Mx
            My, //nullptr, // My
            Mz //nullptr  // Mz
        );

        mx = mean(Mx, Nx*Ny);
        my = mean(My, Nx*Ny);
        mz = mean(Mz, Nx*Ny);

        normalize(mx, my, mz);

        fout << Bext << " " << mx << " " << my << " " << mz<<"\n"; 
        // cout << Bext << " " << mx << " " << my << " " << mz<<"\n";  
        cout<< max_counter-- <<"\n"; 
    }

    //  3rd part of the Hysteresis loop    
    for (int l=-(Bmax/dB); l<= (Bmax/dB); l++){  // for (double Bext = -Bmax; Bext <= Bmax; Bext += dB){
        Bext = l*dB;
        llg_solver.solve(
            // Physics
            50.0,   // t_switch
            Bext,    // Bext
            Aexch,    // Aexch
            D,    // D
            alpha,   // alpha
            gamma_gyro,    // gamma_gyro

            // grid/time
            Nt,      // Nt
            Nx,     // Nx
            Ny,     // Ny
            t,    // t
            dt,   // dt

            // arrays
            Mx, //nullptr, // Mx
            My, //nullptr, // My
            Mz //nullptr  // Mz
        );

        mx = mean(Mx, Nx*Ny);
        my = mean(My, Nx*Ny);
        mz = mean(Mz, Nx*Ny);

        normalize(mx, my, mz);

        fout << Bext << " " << mx << " " << my << " " << mz<<"\n"; 
        // cout << Bext << " " << mx << " " << my << " " << mz<<"\n";  
        cout<< max_counter-- <<"\n";
    }

 
   
    fout.close(); // Close the file
    

    // delete heap memory allocation
    delete[] Mx;
    delete[] My;
    delete[] Mz;

    return 0;
}

   