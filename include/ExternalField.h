/*
===============================================================================
ExternalField Class
-------------------------------------------------------------------------------
Provides different time-dependent external magnetic field configurations
for use in LLG simulations.

All fields are applied along the z-direction only:
    H = (0, 0, Hz)

Sign convention:
    The field is applied as Hz = -B(t)
    (negative sign consistent with chosen energy/LLG convention).

Implemented field types:

1) Static DC Field
   - Time independent
   - Hz = -B0

2) AC Field
   - Sinusoidal time dependence
   - Hz = -B0 sin(omega t)

3) Step Field
   - Constant field before time t0
   - Switched off after t0
   - Hz = -B0  for t < t0
   - Hz = 0    for t >= t0

All functions:
    - Take time t as input
    - Return field components by reference
    - Do not store internal state
    - Designed for lightweight inline evaluation inside time integrators

Units:
    - B0 must be consistent with the chosen gamma and time units
    - omega in angular frequency units
===============================================================================
*/

#ifndef EXTERNALFIELD_H
#define EXTERNALFIELD_H

#include <cmath>

class ExternalField{
    public:
        ExternalField() = default;

        // ---------- Static DC Field ---------
        inline void staticB(double t, double B0, double &Hx, double &Hy, double &Hz){
            Hx = 0.0;
            Hy = 0.0;
            Hz = -B0;
        }

        // ---------- AC Field ----------
        inline void acB(double t, double B0, double omega,  
            double &Hx, double &Hy, double &Hz){
            Hx = 0.0;
            Hy = 0.0;
            Hz = -B0*std::sin(omega*t);
        }

        // ---------- Step Field ----------
        inline void stepB(double t, double t0, double B0,  
            double &Hx, double &Hy, double &Hz){
            Hx = 0.0;
            Hy = 0.0;
            Hz = (t<t0) ? -B0 : 0;
        }
};

#endif


// #include "ExternalField.h"
// #include <iostream>

// int main() {
//     double t = 0.0;        // time
//     double B0 = 3.0;       // field magnitude

//     double h1x, h1y, h1z;

//     ExternalField Bext;
//     Bext.staticB(t, B0, h1x, h1y, h1z);

//     std::cout << h1z << std::endl;    
//     return 0;
// }