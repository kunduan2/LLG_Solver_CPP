#include "utils.h"
#include "ExternalField.h"
#include "ExchangeField.h"
#include "DelmLlg.h"
#include "LlgSolver.h"

struct MaterialParameters{
    double Bext;
    double Aexch;
    double randD;
    double alpha;
    double gamma_gyro;
};

// struct GridTimeParams{
//     int Nt;
//     int Nx;
//     int Ny;
//     double t;
//     double dt;
//     double t_switch;
// };

// struct MagnetizationVectors{
//     double* Mx;    
//     double* My;     
//     double* Mz;   
// };

// std::mt19937 gen{42};
// std::normal_distribution<double> gauss{0.0, 1.0};

void LlgSolver::solve(
    const MaterialParameters& matparams,
    
    // grid/time
    int Nt,
    int Nx,
    int Ny,
    double t,
    double dt,
    
    double* Mx,      
    double* My,     
    double* Mz  

){                
        double Hex, Hey, Hez;
        double Hext1x, Hext1y, Hext1z, Hext2x, Hext2y, Hext2z, Hext3x, Hext3y, Hext3z, Hext4x, Hext4y, Hext4z;
        double k1x, k1y, k1z, k2x, k2y, k2z, k3x, k3y, k3z, k4x, k4y, k4z;
        double Htotx, Htoty, Htotz;  

        // Define objects
        ExchangeField exchange_field(Nx, Ny);
        ExternalField external_field;
        DelmLlg dmdt;

        for (int tn = 1; tn < Nt; tn++) {
            for (int i = 0; i < Ny - 1; i++) {
                for (int j = 0; j < Nx - 1; j++) {

                    // Exchange-field
                    exchange_field.calculate(
                        i, j, 
                        Mx, My, Mz,
                        Hex, Hey, Hez
                    );

                    // 2D to 1D array-index conversion
                    int idx = convert_idx_2dto1D(i, j, Nx, Ny);
                    
                    double mx = Mx[idx];
                    double my = My[idx];
                    double mz = Mz[idx];
                    
                    /*--------------- add external-field and calculate delm using RK4 ------------------*/
                    // step: RK4-1
                    external_field.staticB(
                        t, matparams.Bext,
                        Hext1x, Hext1y, Hext1z // output variables
                    );

                    Htotx = Hext1x + matparams.Aexch*Hex; // total field (witout the random field)
                    Htoty = Hext1y + matparams.Aexch*Hey;
                    Htotz = Hext1z + matparams.Aexch*Hez;

                    dmdt.calculate(
                        mx, my, mz, 
                        Htotx, Htoty, Htotz,
                        matparams.randD, dt, matparams.alpha, matparams.gamma_gyro,
                        gen, gauss, 
                        k1x, k1y, k1z  // output variables
                    );
                    
                    // step: RK4-2
                    external_field.staticB(
                        t + 0.5*dt, matparams.Bext,
                        Hext2x, Hext2y, Hext2z // output variables
                    );

                    Htotx = Hext2x + matparams.Aexch*Hex; // total field (witout the random field)
                    Htoty = Hext2y + matparams.Aexch*Hey;
                    Htotz = Hext2z + matparams.Aexch*Hez;

                    dmdt.calculate(
                        mx + 0.5*dt*k1x, my + 0.5*dt*k1y, mz + 0.5*dt*k1z,
                        Htotx, Htoty, Htotz,
                        matparams.randD, dt, matparams.alpha, matparams.gamma_gyro, 
                        gen, gauss, 
                        k2x, k2y, k2z // output variables
                    );
                    
                    // step: RK4-3
                    external_field.staticB(
                        t + 0.5*dt, matparams.Bext, 
                        Hext3x, Hext3y, Hext3z // output variables
                    );

                    Htotx = Hext3x + matparams.Aexch*Hex; // total field (witout the random field)
                    Htoty = Hext3y + matparams.Aexch*Hey;
                    Htotz = Hext3z + matparams.Aexch*Hez;

                    dmdt.calculate(
                        mx + 0.5*dt*k2x, my + 0.5*dt*k2y, mz + 0.5*dt*k2z,
                        Htotx, Htoty, Htotz, 
                        matparams.randD, dt, matparams.alpha, matparams.gamma_gyro,
                        gen, gauss,
                        k3x, k3y, k3z // output variables
                    );

                    
                    // step: RK4-4
                    external_field.staticB(
                        t + dt, matparams.Bext, 
                        Hext4x, Hext4y, Hext4z // output variables
                    );

                    Htotx = Hext4x + matparams.Aexch*Hex; // total field (witout the random field)
                    Htoty = Hext4y + matparams.Aexch*Hey;
                    Htotz = Hext4z + matparams.Aexch*Hez;    
                                        
                    dmdt.calculate(
                        mx + dt*k3x, my + dt*k3y, mz + dt*k3z,
                        Htotx, Htoty, Htotz,
                        matparams.randD, dt, matparams.alpha, matparams.gamma_gyro,
                        gen, gauss, 
                        k4x, k4y, k4z
                    ); 
                    
                    // m(t+1) = m(t) + f(m)dt   =  m(t) + (dmdt)dt                  
                    mx += (dt/6.0)*(k1x + 2*k2x + 2*k3x + k4x);
                    my += (dt/6.0)*(k1y + 2*k2y + 2*k3y + k4y);
                    mz += (dt/6.0)*(k1z + 2*k2z + 2*k3z + k4z);

                    normalize(mx, my, mz);

                    // update mi's
                    Mx[idx] = mx;
                    My[idx] = my;
                    Mz[idx] = mz;

                } // end: j-loop

            } // end: i-loop


            t += dt;

        } // end: t-loop
        
    }





















###############################################################################################################



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
//         ExchangeField exchange_field(Nx, Ny);
//         ExternalField external_field;
//         DelmLlg dmdt;

//         for (int tn = 1; tn < Nt; tn++) {
//             for (int i = 0; i < Ny - 1; i++) {
//                 for (int j = 0; j < Nx - 1; j++) {

//                     // Exchange-field
//                     exchange_field.calculate(
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
//                     external_field.stepB(
//                         t, t_switch, matparams.Bext,
//                         Hext1x, Hext1y, Hext1z // output variables
//                     );

//                     Htotx = Hext1x + matparams.Aexch*Hex; // total field (witout the random field)
//                     Htoty = Hext1y + matparams.Aexch*Hey;
//                     Htotz = Hext1z + matparams.Aexch*Hez;

//                     dmdt.calculate(
//                         mx, my, mz, 
//                         Htotx, Htoty, Htotz,
//                         D, dt, matparams.alpha, matparams.gamma_gyro,
//                         gen, gauss, 
//                         k1x, k1y, k1z  // output variables
//                     );
                    
//                     // step: RK4-2
//                     external_field.stepB(
//                         t + 0.5*dt,  t_switch, matparams.Bext,
//                         Hext2x, Hext2y, Hext2z // output variables
//                     );

//                     Htotx = Hext2x + matparams.Aexch*Hex; // total field (witout the random field)
//                     Htoty = Hext2y + matparams.Aexch*Hey;
//                     Htotz = Hext2z + matparams.Aexch*Hez;

//                     dmdt.calculate(
//                         mx + 0.5*dt*k1x, my + 0.5*dt*k1y, mz + 0.5*dt*k1z,
//                         Htotx, Htoty, Htotz,
//                         D, dt, matparams.alpha, matparams.gamma_gyro, 
//                         gen, gauss, 
//                         k2x, k2y, k2z // output variables
//                     );
                    
//                     // step: RK4-3
//                     external_field.stepB(
//                         t + 0.5*dt,  t_switch, matparams.Bext, 
//                         Hext3x, Hext3y, Hext3z // output variables
//                     );

//                     Htotx = Hext3x + matparams.Aexch*Hex; // total field (witout the random field)
//                     Htoty = Hext3y + matparams.Aexch*Hey;
//                     Htotz = Hext3z + matparams.Aexch*Hez;

//                     dmdt.calculate(
//                         mx + 0.5*dt*k2x, my + 0.5*dt*k2y, mz + 0.5*dt*k2z,
//                         Htotx, Htoty, Htotz, 
//                         D, dt, matparams.alpha, matparams.gamma_gyro,
//                         gen, gauss,
//                         k3x, k3y, k3z // output variables
//                     );

                    
//                     // step: RK4-4
//                     external_field.stepB(
//                         t + dt, t_switch, matparams.Bext, 
//                         Hext4x, Hext4y, Hext4z // output variables
//                     );

//                     Htotx = Hext4x + matparams.Aexch*Hex; // total field (witout the random field)
//                     Htoty = Hext4y + matparams.Aexch*Hey;
//                     Htotz = Hext4z + matparams.Aexch*Hez;    
                                        
//                     dmdt.calculate(
//                         mx + dt*k3x, my + dt*k3y, mz + dt*k3z,
//                         Htotx, Htoty, Htotz,
//                         D, dt, matparams.alpha, matparams.gamma_gyro,
//                         gen, gauss, 
//                         k4x, k4y, k4z
//                     ); 
                    
//                     // m(t+1) = m(t) + f(m)dt   =  m(t) + (dmdt)dt                  
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
