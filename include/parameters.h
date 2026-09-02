#ifndef PARAMETERS_H
#define PARAMETERS_H

struct MaterialParameters {
    double Bext       = 0.0;   // External magnetic field
    double Aexch      = 0.0;   // Exchange stiffness
    double D          = 0.0;   // Noise amplitude
    double alpha      = 0.0;   // Damping coefficient
    double gamma_gyro = 1.0;   // Gyromagnetic ratio
    double Dmi_x      = 0.0;   // DMI vector x-component
    double Dmi_y      = 0.0;   // DMI vector y-component
    double Dmi_z      = 0.0;   // DMI vector z-component
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