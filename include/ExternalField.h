
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



// sample test code

// #include "ExternalField.h"
// #include <iostream>
// #include <cmath>
// using namespace std;

// int main() {
//     double h = 3;
//     double h1x, h1y, h1z;

//     ExternalField Bext;
//     Bext.staticB(h, h1x, h1y, h1z) ;

//     std::cout<<h1z << std::endl;    
//     return 0;
// }
