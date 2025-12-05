#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

/*
At T = 0K, no stochastic field \eta = 0 
*/

// Physical parameters for the Landau-Lifshitz-Gilbert equation
const double gamma_gyro = 1.0;  // Gyromagnetic ratio - determines precession speed
const double alpha = 0.1;       // Gilbert damping parameter 

// External magnetic field components 
const double Hx = 0.0;  
const double Hy = 0.0;  
const double Hz = -0.7;  

// Time parameters for numerical integration
double t = 0.0;        
const double dt = 0.01;      
const int N = 10000;  

// Function to compute time derivatives of magnetization using LLG equation
// Input: current magnetization components (mx, my, mz) - should be unit vector
// Output: time derivatives (dmx, dmy, dmz) - passed by reference
void dmdt(double mx, double my, double mz, double &dmx, double &dmy, double &dmz) {
    // Standard Landau-Lifshitz-Gilbert equation:
    // dm/dt = -γ/(1+α²) [m × H + α(m × (m × H))]
    // First term: precession around H, Second term: damping toward H
    
    double prefactor = -gamma_gyro / (1.0 + alpha*alpha);  // Common prefactor for both terms
    
    // Compute m × H (cross product): precession term
    double cross_x = my*Hz - mz*Hy;  
    double cross_y = mz*Hx - mx*Hz;  
    double cross_z = mx*Hy - my*Hx;  

    // Compute m × (m × H)= H - (m·H)m     
    double m_dot_H = mx*Hx + my*Hy + mz*Hz;  // Dot product m·H
    double double_cross_x = Hx - m_dot_H * mx;  
    double double_cross_y = Hy - m_dot_H * my;  
    double double_cross_z = Hz - m_dot_H * mz; 
    
    // Combine precession and damping terms with prefactor
    dmx = prefactor * (cross_x + alpha * double_cross_x);  
    dmy = prefactor * (cross_y + alpha * double_cross_y);  
    dmz = prefactor * (cross_z + alpha * double_cross_z);  
}

// Function to normalize magnetization vector to unit length, i.e., mx² + my² + mz² = 1 
void normalize_magnetization(double &mx, double &my, double &mz) {   // passed by reference/pointer
    double norm = sqrt(mx*mx + my*my + mz*mz);  // Calculate vector magnitude/Norm
    if (norm > 0) {
        mx /= norm;  
        my /= norm;   
        mz /= norm;  
    }
}

int main() {
    // Initial magnetization vector 
    double mx = 1.0;
    double my = 0.0;
    double mz = 0.0;

    // Make m to be a unit vector
    normalize_magnetization(mx, my, mz);          

    // RK4 (Runge-Kutta 4th order) intermediate slopes    
    double k1x, k1y, k1z;  // Slope at beginning of interval
    double k2x, k2y, k2z;  // Slope at midpoint using k1
    double k3x, k3y, k3z;  // Slope at midpoint using k2  
    double k4x, k4y, k4z;  // Slope at end using k3

    // Open output file to save magnetization
    ofstream fout("magnetization.dat");
    fout << t << " " << mx << " " << my << " " << mz << "\n";  // Write initial state

    // Evolve the system over time
    for (int i = 0; i < N; i++) {
        // RK4 Step 1: Compute k1 (slope at beginning) ---
        dmdt(mx, my, mz, k1x, k1y, k1z);  // Get all k1's as output
        
        // RK4 Step 2: Compute k2 (slope at midpoint using k1) 
        dmdt(mx + 0.5*dt*k1x, my + 0.5*dt*k1y, mz + 0.5*dt*k1z, k2x, k2y, k2z); // Get all k2's as output
        
        // RK4 Step 3: Compute k3 (slope at midpoint using k2) 
        dmdt(mx + 0.5*dt*k2x, my + 0.5*dt*k2y, mz + 0.5*dt*k2z, k3x, k3y, k3z); // Get all k3's as output
        
        // RK4 Step 4: Compute k4 (slope at end using k3) 
        dmdt(mx + dt*k3x, my + dt*k3y, mz + dt*k3z, k4x, k4y, k4z); // Get all k4's as output

        // RK4 Step 5: Combine slopes to update magnetization 
        mx += dt/6.0*(k1x + 2*k2x + 2*k3x + k4x); 
        my += dt/6.0*(k1y + 2*k2y + 2*k3y + k4y); 
        mz += dt/6.0*(k1z + 2*k2z + 2*k3z + k4z); 

        // Normalize magnetization to preserve physical constraint |m| = 1
        normalize_magnetization(mx, my, mz);

        // Advance time
        t += dt;

        // Save current state to file for analysis/plotting
        fout << t << " " << mx << " " << my << " " << mz << "\n";
    }

    // Clean up and close output file
    fout.close();
    cout << "LLG simulation complete. Data saved in 'magnetization.dat'\n";
    
    // Verify that magnetization magnitude is preserved (should be ~1.0)
    double final_norm = sqrt(mx*mx + my*my + mz*mz);
    cout << "Final magnetization magnitude: " << final_norm << " (should be close to 1.0)\n";
    
    return 0;  // Successful program termination
}