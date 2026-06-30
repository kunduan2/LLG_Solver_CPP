#include "utils.h"
#include "ExternalField.h"
#include "ExchangeField.h"
#include "LlgStep.h"
#include "LlgSolver.h"

#include <iostream>
#include <random>
#include <fstream>
using namespace std;



int main() {

    // -------------------------------------------------------------------------
    // Define Material parameters
    // -------------------------------------------------------------------------
    const double Aexch      = 1.0;    // Exchange stiffness
    const double gamma_gyro = 0.10;   // Gyromagnetic ratio
    const double alpha      = 0.2;    // Gilbert damping coefficient
    const double kB         = 1.0;    // Boltzmann constant (reduced units)
    const double T          = 0.0;  // Temperature [K]
    const double Ms         = 1.0;    // Saturation magnetization
    const double V          = 1.0;    // Magnetic moment volume
    const double D          = alpha * kB * T / (gamma_gyro * Ms * V); // Thermal noise amplitude

    // -------------------------------------------------------------------------
    //  Set the struct: MaterialParameters (from above values) and  Simulation parameters
    // -------------------------------------------------------------------------
   
    MaterialParameters matparams{
        1.0,        // Bext       : external field
        Aexch,      // Aexch      : exchange stiffness
        D,          // D          : stochastic noise amplitude
        alpha,      // alpha      : damping
        gamma_gyro  // gamma_gyro : gyromagnetic ratio
    };
    
    SimulationParameters simparams{
        100, 100,   // Nx, Ny  : grid dimensions
        50,       // Nt      : number of time steps
        0.01,   // dt      : time step size
    };

    // -------------------------------------------------------------------------
    // RNG
    // -------------------------------------------------------------------------
    /*
    One generator for the whole simulation, declared outside any loop so it's
    constructed once and its state persists across calls (do not reseed).
    If wrapped inside a function, declare it `static` so the same instance is
    reused across calls instead of being reconstructed each time.
    Likewise, any std::normal_distribution used with it should be declared
    `static` in single-threaded code to preserve its internal cache.
    */
    std::mt19937 gen(42);

    // -------------------------------------------------------------------------
    // Make magnetization grid and initialize
    // -------------------------------------------------------------------------
    double* Mx = new double[simparams.Nx * simparams.Ny];
    double* My = new double[simparams.Nx * simparams.Ny];
    double* Mz = new double[simparams.Nx * simparams.Ny];

    // Open a file
    ofstream fout("test.dat");
    if (!fout.is_open()) {
        cerr << "Error opening file test.dat!" << endl;
        return 1; // or handle error appropriately
    }

    // initiatlize magnetization: random state
    for (int i=0; i<simparams.Ny; i++){
        for (int j=0; j<simparams.Nx; j++){
            int ind_ij = convert_idx_2dto1D(i, j, simparams.Nx, simparams.Ny);            
            Mx[ind_ij]= rand_uniform(gen, -1.0, 1.0);
            My[ind_ij]= rand_uniform(gen, -1.0, 1.0);
            Mz[ind_ij]= rand_uniform(gen, -1.0, 1.0);

            normalize(Mx[ind_ij], My[ind_ij], Mz[ind_ij]);
            // cout << Mz[ind_ij] << " ";
        }
    }

    // save mean vals of initiatlize magnetization components
    fout << 0 << " " 
        << mean(Mx, simparams.Nx*simparams.Ny) 
        << " " << mean(My, simparams.Nx*simparams.Ny)
        << " " << mean(Mz, simparams.Nx*simparams.Ny)
        <<"\n";   
    
    // -------------------------------------------------------------------------
    // LLG_Solver
    // -------------------------------------------------------------------------
    LlgSolver llg_solver;
   
    // -------------------------------------------------------------------------
    // Hysteresis-loop
    // -------------------------------------------------------------------------
    double Bmax = 10.2;
    double dB = .50;
    int max_counter = 5*(Bmax/dB); // To track the loop  

    // define variables to be used in all for loops
    double mx = 0.0;
    double my = 0.0;
    double mz = 0.0;
    double Bext = 0.0;

    // 1st part of the Hysteresis loop
    for (int l=0; l<= (Bmax/dB); l++){
        Bext = l*dB; 
        matparams.Bext = Bext;   // update before each solve         
        llg_solver.solve(matparams, simparams, gen, Mx, My, Mz);

        mx = mean(Mx, simparams.Nx*simparams.Ny);
        my = mean(My, simparams.Nx*simparams.Ny);
        mz = mean(Mz, simparams.Nx*simparams.Ny);

        normalize(mx, my, mz);

        fout << Bext << " " << mx << " " << my << " " << mz<<"\n"; 
        // cout << Bext << " " << mx << " " << my << " " << mz<<"\n";  
        cout<< max_counter-- <<"\n";
    }

     // 2nd part of the Hysteresis loop
    for (int l=(Bmax/dB); l>= -(Bmax/dB); l--){  // for (double Bext = Bmax; Bext >= -Bmax; Bext -= dB){
        Bext = l*dB;  
        matparams.Bext = Bext;   // update before each solve         
        llg_solver.solve(matparams, simparams, gen, Mx, My, Mz);

        mx = mean(Mx, simparams.Nx*simparams.Ny);
        my = mean(My, simparams.Nx*simparams.Ny);
        mz = mean(Mz, simparams.Nx*simparams.Ny);

        normalize(mx, my, mz);

        fout << Bext << " " << mx << " " << my << " " << mz<<"\n"; 
        // cout << Bext << " " << mx << " " << my << " " << mz<<"\n";  
        cout<< max_counter-- <<"\n"; 
    }

    //  3rd part of the Hysteresis loop    
    for (int l=-(Bmax/dB); l<= (Bmax/dB); l++){  // for (double Bext = -Bmax; Bext <= Bmax; Bext += dB){
        Bext = l*dB;
        matparams.Bext = Bext;   // update before each solve         
        llg_solver.solve(matparams, simparams, gen, Mx, My, Mz);

        mx = mean(Mx, simparams.Nx*simparams.Ny);
        my = mean(My, simparams.Nx*simparams.Ny);
        mz = mean(Mz, simparams.Nx*simparams.Ny);

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

   