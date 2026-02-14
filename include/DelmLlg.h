#ifndef DELMLLG_H
#define DELMLLG_H

#include <random>

class DelmLlg{
public:
    DelmLlg() = default;

    void calculate(
            // mi's
        double mx, 
        double my, 
        double mz,

        // total external fields
        double Htx, 
        double Hty, 
        double Htz,

        //parameters 
        double D, 
        double dt, 
        double alpha, 
        double gamma_gyro,

        // RNG
        std::mt19937 &gen,
        std::normal_distribution<double> &gauss,

        // output variables
        double &dmx, double &dmy, double &dmz
    );

};


#endif



// sample code to test

// #include "DelmLlg.h"
// #include <iostream>
// #include <random>

// int main() {
//     // ---------- Magnetization (unit vector) ----------
//     double mx = 0.0;
//     double my = 0.0;
//     double mz = 1.0;

//     // ---------- Effective field ----------
//     double Hx = 0.0;
//     double Hy = 0.0;
//     double Hz = 1.0;

//     // ---------- LLG parameters ----------
//     double D = 0.01;
//     double dt = 1e-3;
//     double alpha = 0.1;
//     double gamma = 1.0;

//     // ---------- RNG ----------
//     std::mt19937 gen(42);  // reproducible
//     std::normal_distribution<double> gauss(0.0, 1.0);

//     // ---------- Output ----------
//     double dmx, dmy, dmz;

//     // ---------- LLG object ----------
//     DelmLlg llg;

//     // ---------- Single evaluation ----------
//     llg.calculate(mx, my, mz,
//                   Hx, Hy, Hz,
//                   D, dt, alpha, gamma,
//                   gen, gauss,
//                   dmx, dmy, dmz);

//     std::cout << "dm = ("
//               << dmx << ", "
//               << dmy << ", "
//               << dmz << ")\n";

//     return 0;
// }
