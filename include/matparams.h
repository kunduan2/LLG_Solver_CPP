#ifndef MATPARAMS_H
#define MATPARAMS_H

struct MaterialParameters {
    double Bext;        // External magnetic field
    double Aexch;       // Exchange stiffness
    double D;           // Noise amplitude
    double alpha;       // Damping coefficient
    double gamma_gyro;  // Gyromagnetic ratio
    double D_x;       // DMI vector x-component
    double D_y;       // DMI vector y-component
    double D_z;       // DMI vector z-component
};

struct SimulationParameters {
    int Nx, Ny;         // Grid dimensions
    int Nt;             // Number of time steps
    double dt;          // Time step
    // double t_start;     // Initial time
    
    // // Temperature parameters
    // double T;           // Temperature
    // double kB;          // Boltzmann constant
    // double Ms;          // Saturation magnetization
    // double V;           // Volume
};

#endif // MATPARAMS_H