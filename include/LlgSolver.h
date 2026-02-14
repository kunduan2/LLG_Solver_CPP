#ifndef LLGSOLVER_H
#define LLGSOLVER_H

#include <random>
// #include <fstream>

struct MaterialParameters{
    double Bext;
    double Aexch;
    double randD;
    double alpha;
    double gamma_gyro;
};

// Solves LLG equation for all spins over time using Heun's/RK4 method
// Simulates hysteresis by switching external field at t_switch
class LlgSolver{
private:
    std::mt19937 gen{42};
    std::normal_distribution<double> gauss{0.0, 1.0};
public:
void solve(
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
);
};


#endif