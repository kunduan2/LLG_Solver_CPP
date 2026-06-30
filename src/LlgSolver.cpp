#include <iostream>
#include "utils.h"
#include "matparams.h"

// Physical field calculations
#include "ExchangeField.h"      // Exchange interaction field (ferromagnetic coupling)
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
    double Hex, Hey, Hez;
    double Hext1x, Hext1y, Hext1z;  // Hext2x, Hext2y, Hext2z, Hext3x, Hext3y, Hext3z, Hext4x, Hext4y, Hext4z;
    double k1x, k1y, k1z, k2x, k2y, k2z, k3x, k3y, k3z, k4x, k4y, k4z;
    double Htotx, Htoty, Htotz;  
    double etax, etay, etaz; // Stochastic/thermal fields

    // Define objects
    ExchangeField hexch(simparams.Nx, simparams.Ny);
    ExternalField hext;  // We use static B_ext
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
                Htotx = matparams.Aexch*Hex + Hext1x + etax; 
                Htoty = matparams.Aexch*Hey + Hext1y + etay;
                Htotz = matparams.Aexch*Hez + Hext1z + etaz;

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












// // std::mt19937 gen{42};
// // std::normal_distribution<double> gauss{0.0, 1.0};

// void LlgSolver::solve(

//         // Physics
//         double t_switch,
//         double matparams.Bext,
//         double matparams.Aexch,
//         double D ,
//         double matparams.alpha,
//         double matparams.gamma_gyro,

//         // grid/time
//         int Nt,
//         int Nx,
//         int Ny,
//         double t,
//         double dt,
        
//         double* Mx,      
//         double* My,     
//         double* Mz      
//     ){                
//         double Hex, Hey, Hez;
//         double Hext1x, Hext1y, Hext1z, Hext2x, Hext2y, Hext2z, Hext3x, Hext3y, Hext3z, Hext4x, Hext4y, Hext4z;
//         double k1x, k1y, k1z, k2x, k2y, k2z, k3x, k3y, k3z, k4x, k4y, k4z;
//         double Htotx, Htoty, Htotz;  

//         // Define objects
//         ExchangeField hexch(Nx, Ny);
//         ExternalField hext;
//         DelmLlg delm;

//         for (int tn = 1; tn < Nt; tn++) {
//             for (int i = 0; i < Ny - 1; i++) {
//                 for (int j = 0; j < Nx - 1; j++) {

//                     // Exchange-field
//                     hexch.calculate(
//                         i, j, 
//                         Mx, My, Mz,
//                         Hex, Hey, Hez
//                     );

//                     // 2D to 1D array-index conversion
//                     int idx = convert_idx_2dto1D(i, j, Nx, Ny);
                    
//                     double mx = Mx[idx];
//                     double my = My[idx];
//                     double mz = Mz[idx];
                    
//                     /*--------------- add external-field and calculate delm using RK4 ------------------*/
//                     // step: RK4-1
//                     hext.stepB(
//                         t, t_switch, matparams.Bext,
//                         Hext1x, Hext1y, Hext1z // output variables
//                     );

//                     Htotx = Hext1x + matparams.Aexch*Hex; // total field (witout the random field)
//                     Htoty = Hext1y + matparams.Aexch*Hey;
//                     Htotz = Hext1z + matparams.Aexch*Hez;

//                     delm.calculate(
//                         mx, my, mz, 
//                         Htotx, Htoty, Htotz,
//                         D, dt, matparams.alpha, matparams.gamma_gyro,
//                         gen, gauss, 
//                         k1x, k1y, k1z  // output variables
//                     );
                    
//                     // step: RK4-2
//                     hext.stepB(
//                         t + 0.5*dt,  t_switch, matparams.Bext,
//                         Hext2x, Hext2y, Hext2z // output variables
//                     );

//                     Htotx = Hext2x + matparams.Aexch*Hex; // total field (witout the random field)
//                     Htoty = Hext2y + matparams.Aexch*Hey;
//                     Htotz = Hext2z + matparams.Aexch*Hez;

//                     delm.calculate(
//                         mx + 0.5*dt*k1x, my + 0.5*dt*k1y, mz + 0.5*dt*k1z,
//                         Htotx, Htoty, Htotz,
//                         D, dt, matparams.alpha, matparams.gamma_gyro, 
//                         gen, gauss, 
//                         k2x, k2y, k2z // output variables
//                     );
                    
//                     // step: RK4-3
//                     hext.stepB(
//                         t + 0.5*dt,  t_switch, matparams.Bext, 
//                         Hext3x, Hext3y, Hext3z // output variables
//                     );

//                     Htotx = Hext3x + matparams.Aexch*Hex; // total field (witout the random field)
//                     Htoty = Hext3y + matparams.Aexch*Hey;
//                     Htotz = Hext3z + matparams.Aexch*Hez;

//                     delm.calculate(
//                         mx + 0.5*dt*k2x, my + 0.5*dt*k2y, mz + 0.5*dt*k2z,
//                         Htotx, Htoty, Htotz, 
//                         D, dt, matparams.alpha, matparams.gamma_gyro,
//                         gen, gauss,
//                         k3x, k3y, k3z // output variables
//                     );

                    
//                     // step: RK4-4
//                     hext.stepB(
//                         t + dt, t_switch, matparams.Bext, 
//                         Hext4x, Hext4y, Hext4z // output variables
//                     );

//                     Htotx = Hext4x + matparams.Aexch*Hex; // total field (witout the random field)
//                     Htoty = Hext4y + matparams.Aexch*Hey;
//                     Htotz = Hext4z + matparams.Aexch*Hez;    
                                        
//                     delm.calculate(
//                         mx + dt*k3x, my + dt*k3y, mz + dt*k3z,
//                         Htotx, Htoty, Htotz,
//                         D, dt, matparams.alpha, matparams.gamma_gyro,
//                         gen, gauss, 
//                         k4x, k4y, k4z
//                     ); 
                    
//                     // m(t+1) = m(t) + f(m)dt   =  m(t) + (delm)dt                  
//                     mx += (dt/6.0)*(k1x + 2*k2x + 2*k3x + k4x);
//                     my += (dt/6.0)*(k1y + 2*k2y + 2*k3y + k4y);
//                     mz += (dt/6.0)*(k1z + 2*k2z + 2*k3z + k4z);

//                     normalize(mx, my, mz);

//                     // update mi's
//                     Mx[idx] = mx;
//                     My[idx] = my;
//                     Mz[idx] = mz;

//                 } // end: j-loop

//             } // end: i-loop


//             t += dt;

//         } // end: t-loop
        
//     }
