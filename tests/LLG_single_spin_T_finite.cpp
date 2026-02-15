#include <iostream>
#include <fstream>
#include <cmath>
#include <random>   // For Gaussian noise
using namespace std;

// Physical parameters for LLG
const double gamma_gyro = 1.0;
const double alpha = 0.1;
const double prefactor = -gamma_gyro / (1.0 + alpha*alpha);

// Time-loop parameters
const double dt = 0.01;
const int N = 10000;

// --------------------------------------
//  External field H(t): Step function: H changes abruptly at t = t_switch
// --------------------------------------

// Step field parameters
const double H_before = -0.7;    // value before switching time
const double H_after  = -0.7;    // value after switching time
const double t_switch = 25.0;    // switching time in simulation units
 

// Temperature / thermal noise
const double kB = 1.0;      // Boltzmann constant (choose units)
const double T = 800.00;      // Temperature
const double Ms = 1.0;     // Saturation magnetization
const double V = 1.0;      // Volume
const double D = alpha*kB*T/(gamma_gyro*Ms*V); // Noise amplitude

// Random number generator for Gaussian noise
std::mt19937 gen(42); // Fixed seed for reproducibility
std::normal_distribution<double> gauss(0.0, 1.0); // Standard normal

// normalize
void normalize(double &mx, double &my, double &mz)
{
    double n = sqrt(mx*mx + my*my + mz*mz);
    mx /= n; my /= n; mz /= n;
}

void H_of_t(double t, double &Hx, double &Hy, double &Hz){
    // Step function: H changes abruptly at t = t_switch
    if (t < t_switch) {
        Hx = 0.0;
        Hy = 0.0;
        Hz = H_before;
    } 
    else {
        Hx = 0.0;
        Hy = 0.0;
        Hz = H_after;
    }
}

// Function to compute time derivatives of magnetization using LLG equation
void dmdt(double mx, double my, double mz,
          double Hx, double Hy, double Hz,
          double &dmx, double &dmy, double &dmz) {
    // Compute thermal field increments
    double eta_x = sqrt(2.0*D*dt)*gauss(gen);
    double eta_y = sqrt(2.0*D*dt)*gauss(gen);
    double eta_z = sqrt(2.0*D*dt)*gauss(gen);

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



int main() {
    // start time 
    double t = 0.0; 

    // magnatization component at t=0.0
    double mx = 1.0, my = 0.0, mz = 0.0;
    normalize(mx, my, mz);

    // --- compute H at stage times (exact RK4 bookkeeping) ---
    double H1x,H1y,H1z;
    double H2x,H2y,H2z;
    double H3x,H3y,H3z;
    double H4x,H4y,H4z;

    // --- compute K's at stage times (exact RK4 bookkeeping) ---
    double k1x,k1y,k1z, k2x,k2y,k2z, k3x,k3y,k3z, k4x,k4y,k4z;    

    ofstream fout("test.dat");
    fout << t << " " << mx << " " << my << " " << mz << "\n";

    for (int i = 0; i < N; i++) {

        // --- compute H at stage times (exact RK4 bookkeeping) ---
        H_of_t(t,               H1x,H1y,H1z);
        H_of_t(t + 0.5*dt,      H2x,H2y,H2z);
        H_of_t(t + 0.5*dt,      H3x,H3y,H3z);
        H_of_t(t + dt,          H4x,H4y,H4z);

        // generate noise ONCE per time step
        double eta_x = sqrt(2.0*D/dt)*gauss(gen);
        double eta_y = sqrt(2.0*D/dt)*gauss(gen);
        double eta_z = sqrt(2.0*D/dt)*gauss(gen);
        
        // RK4 stages (pass stage H into dmdt)
        dmdt(mx, my, mz, 
            H1x,H1y,H1z, k1x,k1y,k1z);

        dmdt(mx + 0.5*dt*k1x, my + 0.5*dt*k1y, mz + 0.5*dt*k1z,
             H2x, H2y, H2z, k2x, k2y, k2z);

        dmdt(mx + 0.5*dt*k2x, my + 0.5*dt*k2y, mz + 0.5*dt*k2z,
             H3x,H3y,H3z, k3x,k3y,k3z);

        dmdt(mx + dt*k3x, my + dt*k3y, mz + dt*k3z,
             H4x,H4y,H4z, k4x,k4y,k4z);

        mx += dt/6.0*(k1x + 2*k2x + 2*k3x + k4x);
        my += dt/6.0*(k1y + 2*k2y + 2*k3y + k4y);
        mz += dt/6.0*(k1z + 2*k2z + 2*k3z + k4z);

        normalize(mx, my, mz);
        t += dt;
        fout << t << " " << mx << " " << my << " " << mz << "\n";
    }

    fout.close();
    // cout << "Stochastic LLG simulation complete. Data saved in 'magnetization_stochastic.dat'\n";
    // double final_norm = sqrt(mx*mx + my*my + mz*mz);
    // cout << "Final magnetization magnitude: " << final_norm << " (should be close to 1.0)\n";
    return 0;
}
