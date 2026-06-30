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
    // Material parameters
    // -------------------------------------------------------------------------
    const double kB        = 1.0;    // Boltzmann constant (reduced units)
    const double T         = 100.0;    // Temperature
    const double Ms        = 1.0;    // Saturation magnetization
    const double V         = 1.0;    // Magnetic moment volume
    const double alpha     = .001;    // Gilbert damping coefficient
    const double gamma_gyro= 1.0;    // Gyromagnetic ratio
    const double Aexch     = 0.250;  // Exchange stiffness
    const double D         = alpha * kB * T / (gamma_gyro * Ms * V); // Noise amplitude

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
        5, 5,   // Nx, Ny  : grid dimensions
        5,     // Nt      : number of time steps
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
    // Initialize magnetization grid
    // -------------------------------------------------------------------------
    double* Mx = new double[simparams.Nx * simparams.Ny];
    double* My = new double[simparams.Nx * simparams.Ny];
    double* Mz = new double[simparams.Nx * simparams.Ny];

    // initiatlize magnetization
    double mx = 1.0;
    double my = 0.0;
    double mz = 1.0;
    
    normalize(mx, my, mz);

    for (int i=0; i<simparams.Ny; i++){
        for (int j=0; j<simparams.Nx; j++){
            int ind_ij = convert_idx_2dto1D(i, j, simparams.Nx, simparams.Ny);            
            Mx[ind_ij]= mx;
            My[ind_ij]= my;
            Mz[ind_ij]= mz;
            // cout << Mz[ind_ij] << " ";
        }
    }

//    print2Darr(Mx, simparams.Nx, simparams.Ny);

    LlgSolver llg_solver;
    llg_solver.solve(matparams, simparams, gen, Mx, My, Mz);
    

    // delete heap memory allocation
    delete[] Mx;
    delete[] My;
    delete[] Mz;

    return 0;
}
