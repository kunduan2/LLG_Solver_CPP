#include <iostream>
#include <fstream>
#include <cmath>
#include <random>   // For Gaussian noise
using namespace std;

// Physical parameters for the Landau-Lifshitz-Gilbert equation
const double gamma_gyro = 1.0;  // Gyromagnetic ratio
const double alpha = 0.1;       // Gilbert damping parameter 

// External magnetic field components 
const double Hx = 0.0;  
const double Hy = 0.0;  
const double Hz = -0.7;  

// Temperature / thermal noise
const double kB = 1.0;      // Boltzmann constant (choose units)
const double T = 300;      // Temperature
const double Ms = 1.0;      // Saturation magnetization
const double V = 1.0;       // Volume
const double D = alpha*kB*T/(gamma_gyro*Ms*V); // Noise amplitude

// Time parameters for numerical integration
double t = 0.0;        
const double dt = 0.01;      
const int N = 10000;  

// Random number generator for Gaussian noise
std::mt19937 gen(42); // Fixed seed for reproducibility
std::normal_distribution<double> gauss(0.0, 1.0); // Standard normal

// Function to compute time derivatives of magnetization using LLG equation
void dmdt(double mx, double my, double mz, double &dmx, double &dmy, double &dmz) {
    // Compute thermal field increments
    double eta_x = sqrt(2.0*D*dt)*gauss(gen);
    double eta_y = sqrt(2.0*D*dt)*gauss(gen);
    double eta_z = sqrt(2.0*D*dt)*gauss(gen);

    double prefactor = -gamma_gyro / (1.0 + alpha*alpha);

    // Total field including stochastic term
    double Htot_x = Hx + eta_x;
    double Htot_y = Hy + eta_y;
    double Htot_z = Hz + eta_z;

    // Compute m × Htot (precession term)
    double cross_x = my*Htot_z - mz*Htot_y;  
    double cross_y = mz*Htot_x - mx*Htot_z;  
    double cross_z = mx*Htot_y - my*Htot_x;  

    // Compute m × (m × Htot) (damping term)
    double m_dot_H = mx*Htot_x + my*Htot_y + mz*Htot_z;
    double double_cross_x = Htot_x - m_dot_H * mx;
    double double_cross_y = Htot_y - m_dot_H * my;
    double double_cross_z = Htot_z - m_dot_H * mz;

    // Combine precession and damping
    dmx = prefactor * (cross_x + alpha * double_cross_x);
    dmy = prefactor * (cross_y + alpha * double_cross_y);
    dmz = prefactor * (cross_z + alpha * double_cross_z);
}

// Function to normalize magnetization vector to unit length
void normalize_magnetization(double &mx, double &my, double &mz) {
    double norm = sqrt(mx*mx + my*my + mz*mz);
    if (norm > 0) {
        mx /= norm;
        my /= norm;
        mz /= norm;
    }
}

int main() {
    // initial conditions of mi's
    double mx = 0.0;
    double my = 0.0;
    double mz = 0.0;
    normalize_magnetization(mx, my, mz);

    double k1x, k1y, k1z;
    double k2x, k2y, k2z;
    double k3x, k3y, k3z;
    double k4x, k4y, k4z;

    ofstream fout("magnetization.dat");
    fout << t << " " << mx << " " << my << " " << mz << "\n";

    for (int i = 0; i < N; i++) {
        dmdt(mx, my, mz, k1x, k1y, k1z);
        dmdt(mx + 0.5*dt*k1x, my + 0.5*dt*k1y, mz + 0.5*dt*k1z, k2x, k2y, k2z);
        dmdt(mx + 0.5*dt*k2x, my + 0.5*dt*k2y, mz + 0.5*dt*k2z, k3x, k3y, k3z);
        dmdt(mx + dt*k3x, my + dt*k3y, mz + dt*k3z, k4x, k4y, k4z);

        mx += dt/6.0*(k1x + 2*k2x + 2*k3x + k4x);
        my += dt/6.0*(k1y + 2*k2y + 2*k3y + k4y);
        mz += dt/6.0*(k1z + 2*k2z + 2*k3z + k4z);

        normalize_magnetization(mx, my, mz);
        t += dt;
        fout << t << " " << mx << " " << my << " " << mz << "\n";
    }

    fout.close();
    cout << "Stochastic LLG simulation complete. Data saved in 'magnetization_stochastic.dat'\n";
    double final_norm = sqrt(mx*mx + my*my + mz*mz);
    cout << "Final magnetization magnitude: " << final_norm << " (should be close to 1.0)\n";
    return 0;
}
