#include <iostream>
#include "utils.h"
#include "matparams.h"

// Physical field calculations
#include "ExchangeField.h"      // Exchange interaction field (ferromagnetic coupling)
#include "DMIField.h"      // DMI field
#include "ExternalField.h"      // Applied external magnetic field
#include "StochasticField.h"    // Thermal noise: stochastic field with Gaussian distribution

// Dynamics solver
#include "LlgStep.h"            // Landau-Lifshitz-Gilbert equation formulation
#include "LlgSolver.h"          // LLG equation solver: time integration
// ============================================================================

void LlgSolver::solve(
    const MaterialParameters& matparams,
    const SimulationParameters& simparams, 
    std::mt19937& gen,   
    double* Mx,      
    double* My,     
    double* Mz
){        
    double Hex, Hey, Hez;  // Exchange
    double Hdmi_x, Hdmi_y, Hdmi_z; // DMI
    double Hext1x, Hext1y, Hext1z;  // external
    double k1x, k1y, k1z, k2x, k2y, k2z, k3x, k3y, k3z, k4x, k4y, k4z; // RK-4 steps
    double Htotx, Htoty, Htotz;  // Total 
    double etax, etay, etaz; // Stochastic/thermal fields

    // Define objects
    ExchangeField hexch(simparams.Nx, simparams.Ny);
    ExternalField hext;  // We use static B_ext
    DMIField hdmi(simparams.Nx, simparams.Ny, matparams.Dmi_x, matparams.Dmi_y, matparams.Dmi_z);;
    // StochasticField eta; // No-Need to make instances for static method
    LlgStep delm;

    // // Initialize RNG engine
    // std::mt19937 gen(42);  // seed for reproducibility

    // int itr = 0; // Not used anywhere. just to check code

    double t = 0.0;
    for (int tn = 0; tn < simparams.Nt; tn++) {
        for (int i = 0; i < simparams.Ny; i++) {
            for (int j = 0; j < simparams.Nx ; j++) {
                
                // Exchange-field
                hexch.calculate(
                    i, j, 
                    Mx, My, Mz,
                    Hex, Hey, Hez // output
                );

                // DMI-field
                hdmi.calculate(
                    i, j, 
                    Mx, My, Mz,
                    Hdmi_x, Hdmi_y, Hdmi_z // output
                );

                // stochastic-field
                StochasticField::generate(
                    matparams.D, simparams.dt, gen, 
                    etax, etay, etaz // output
                );
                
                // fixed External-field (loop-independent). t=1 or t=t does not change anything. Check the defn.
                hext.staticB(
                    matparams.Bext,
                    Hext1x, Hext1y, Hext1z // output variables
                );
                
                // total field
                Htotx = matparams.Aexch*Hex + Hdmi_x + Hext1x + etax; 
                Htoty = matparams.Aexch*Hey + Hdmi_y + Hext1y + etay;
                Htotz = matparams.Aexch*Hez + Hdmi_z + Hext1z + etaz;

                // 2D to 1D array-index conversion
                int idx = convert_idx_2dto1D(i, j, simparams.Nx, simparams.Ny);
                
                double mx = Mx[idx];
                double my = My[idx];
                double mz = Mz[idx];
                
                /*--------------- add external-field and calculate delm using RK4 ------------------*/
                // step: RK4-1
                delm.calculate(
                    mx, my, mz, 
                    Htotx, Htoty, Htotz,
                    matparams.alpha, matparams.gamma_gyro,
                    k1x, k1y, k1z  // output variables
                );
                
                // step: RK4-2
                delm.calculate(
                    mx + 0.5*simparams.dt*k1x, my + 0.5*simparams.dt*k1y, mz + 0.5*simparams.dt*k1z,
                    Htotx, Htoty, Htotz,
                    matparams.alpha, matparams.gamma_gyro,
                    k2x, k2y, k2z // output variables
                );
                
                // step: RK4-3 
                delm.calculate(
                    mx + 0.5*simparams.dt*k2x, my + 0.5*simparams.dt*k2y, mz + 0.5*simparams.dt*k2z,
                    Htotx, Htoty, Htotz, 
                    matparams.alpha, matparams.gamma_gyro,
                    k3x, k3y, k3z // output variables
                );

                
                // step: RK4-4                                    
                delm.calculate(
                    mx + simparams.dt*k3x, my + simparams.dt*k3y, mz + simparams.dt*k3z,
                    Htotx, Htoty, Htotz,
                    matparams.alpha, matparams.gamma_gyro,
                    k4x, k4y, k4z
                ); 
                
                // m(t+1) = m(t) + f(m)dt   =  m(t) + (delm)dt                  
                mx += (simparams.dt/6.0)*(k1x + 2*k2x + 2*k3x + k4x);
                my += (simparams.dt/6.0)*(k1y + 2*k2y + 2*k3y + k4y);
                mz += (simparams.dt/6.0)*(k1z + 2*k2z + 2*k3z + k4z);

                normalize(mx, my, mz);

                // update mi's
                Mx[idx] = mx;
                My[idx] = my;
                Mz[idx] = mz;

                // itr += 1; // Not used anywhere. just to check code

                // std::cout<< itr << " " ; // Not used anywhere. just to check code

            } // end: j-loop

        } // end: i-loop

        t += simparams.dt;

    } // end: t-loop
    // std::cout<< "\n\nt = "<< t << "\n";    
}

