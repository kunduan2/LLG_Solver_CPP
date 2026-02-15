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
const int N = 10;



// Temperature / thermal noise
const double kB = 1.0;      // Boltzmann constant (choose units)
const double T = 00.00;      // Temperature
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

// --------------------------------------
//  External field H(t)
// --------------------------------------
const double H_before = 0.7;    // value before switching time
const double H_after  = 0.7;    // value after switching time
const double t_switch = 50.0;    // switching time in simulation units
 
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
// --------------------------------------

// Function to compute time derivatives of magnetization using LLG equation
void delm(double mx, double my, double mz, double Hx, double Hy, double Hz,
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

// void arr3d(double* A, int Nx = 5, int Ny=3){    
//     for (int i=0; i<Ny; i++){
//         for (int j=0; j<Nx; j++){
//             A[i*Nx + j] = i*Nx + j;
//             cout << A[i*Nx + j] << " ";
//         }       
        
//     }
// }

inline int convert_idx_2dto1D (int i, int j, int Nx, int Ny){
    return i*Nx + j ; // Note: Nx is the col numbers and Ny is the row numbers
}

inline int indices_conversion_3dto1D (int i, int j, int k, int Nx, int Ny, int Nz){
    return k*Nx*Ny + i*Nx + j ; // Note: Nx is the col numbers and Ny is the row numbers
}

inline int pbc(int i, int N) {
    return (i + N) % N;
}

void exchange_field (int i, int j, int Nx, int Ny, double* Mx, double* My, double* Mz,
    double &Hexch_x, double &Hexch_y, double &Hexch_z){

    /* --------------------------------------------------------------------------------------------------------------------*/
    /*If no PBC is used, uncomment the following version and comment the PBC consditions next section */
    if (i==0 || j==0 || i==Nx-1 || j==Ny-1){
        Hexch_x =0; Hexch_y =0; Hexch_z=0;
        return;
    }
    int ind_ij   = convert_idx_2dto1D (i, j, Nx, Ny);
    int ind_ip1j = convert_idx_2dto1D (i + 1, j, Nx, Ny);
    int ind_im1j = convert_idx_2dto1D (i - 1, j, Nx, Ny);
    int ind_ijp1 = convert_idx_2dto1D (i, j + 1, Nx, Ny);
    int ind_ijm1 = convert_idx_2dto1D (i, j - 1, Nx, Ny);
    /* --------------------------------------------------------------------------------------------------------------------*/

    // Periodic Boundary conditions (PBC) (OpenMP friendly)
    // int ind_ip1j = convert_idx_2dto1D(ip, j,  Nx, Ny);
    // int ind_im1j = convert_idx_2dto1D(im, j,  Nx, Ny);
    // int ind_ijp1 = convert_idx_2dto1D(i,  jp, Nx, Ny);
    // int ind_ijm1 = convert_idx_2dto1D(i,  jm, Nx, Ny);
 

    // int ind_ij   = convert_idx_2dto1D (i, j, Nx, Ny);
    // int ind_ip1j = convert_idx_2dto1D (i + 1, j, Nx, Ny);
    // int ind_im1j = convert_idx_2dto1D (i - 1, j, Nx, Ny);
    // int ind_ijp1 = convert_idx_2dto1D (i, j + 1, Nx, Ny);
    // int ind_ijm1 = convert_idx_2dto1D (i, j - 1, Nx, Ny);
    /* --------------------------------------------------------------------------------------------------------------------*/

    Hexch_x = Mx[ind_ip1j] + Mx[ind_im1j] + Mx[ind_ijp1] + Mx[ind_ijm1] - 4*Mx[ind_ij];
    Hexch_y = My[ind_ip1j] + My[ind_im1j] + My[ind_ijp1] + My[ind_ijm1] - 4*My[ind_ij];
    Hexch_z = Mz[ind_ip1j] + Mz[ind_im1j] + Mz[ind_ijp1] + Mz[ind_ijm1] - 4*Mz[ind_ij];
}


int main() {
    // start time 
    double t = 0.0; 

    // define constants (will be used inside loop)
    double mx = 0.0, my = 0.0, mz = 0.0;
    int idx;

    // --- compute H at stage times (exact RK4 bookkeeping) ---
    double H1x,H1y,H1z;
    double H2x,H2y,H2z;
    double H3x,H3y,H3z;
    double H4x,H4y,H4z;

    // Exchange 
    double Hexch_x, Hexch_y, Hexch_z;
    double Aexch = 1.0; 

    // Total effective field components (RK4 stages)
    double Htot_RK1_x, Htot_RK1_y, Htot_RK1_z;
    double Htot_RK2_x, Htot_RK2_y, Htot_RK2_z;
    double Htot_RK3_x, Htot_RK3_y, Htot_RK3_z;
    double Htot_RK4_x, Htot_RK4_y, Htot_RK4_z;


    // --- compute K's at stage times (exact RK4 bookkeeping) ---
    double k1x,k1y,k1z, k2x,k2y,k2z, k3x,k3y,k3z, k4x,k4y,k4z;   
    int Nx = 5, Ny=3;
    
    // make a grid for a distribution of magnetic moments in 2D
    double Mx[Nx*Ny];
    double My[Nx*Ny];
    double Mz[Nx*Ny];

    // initiatlize magnetization
    for (int i=0; i<Ny; i++){
        for (int j=0; j<Nx; j++){
            int ind_ij = convert_idx_2dto1D(i, j, Nx, Ny);            
            Mx[ind_ij]= 1;
            My[ind_ij]= 0.01;
            Mz[ind_ij]= 0;
            // cout << Mz[ind_ij] << " ";
        }
    }

    
    for (int tn=0; tn<1000; tn++ ){
        for (int i=0; i<Ny; i++){
            for (int j=0; j<Nx; j++){

            // compute field from Exchange interaction
            exchange_field(i, j, Nx, Ny, Mx, My, Mz,
                        Hexch_x, Hexch_y, Hexch_z); // update the variables using pointers
                
            // --- compute H at stage times (exact RK4 bookkeeping) ---
            H_of_t(t,               H1x, H1y, H1z); // update the variables: H1x, H1y, H1z
            H_of_t(t + 0.5*dt,      H2x, H2y, H2z); 
            H_of_t(t + 0.5*dt,      H3x, H3y, H3z);
            H_of_t(t + dt,          H4x, H4y, H4z);
            
            //-------- Add all the other fields except noise --------------------------------------------
            // X components
            Htot_RK1_x = H1x + Hexch_x;
            Htot_RK2_x = H2x + Hexch_x;
            Htot_RK3_x = H3x + Hexch_x;
            Htot_RK4_x = H4x + Hexch_x;

            // Y components
            Htot_RK1_y = H1y + Hexch_y;
            Htot_RK2_y = H2y + Hexch_y;
            Htot_RK3_y = H3y + Hexch_y;
            Htot_RK4_y = H4y + Hexch_y;

            // Z components
            Htot_RK1_z = H1z + Hexch_z;
            Htot_RK2_z = H2z + Hexch_z;
            Htot_RK3_z = H3z + Hexch_z;
            Htot_RK4_z = H4z + Hexch_z;

            // Mij
            idx =  convert_idx_2dto1D(i, j, Nx, Ny);
            mx = Mx[idx];
            my = My[idx];
            mz = Mz[idx];

            
            // ---------- RK4 stage 1 ----------
            delm(mx, my, mz,
                Htot_RK1_x, Htot_RK1_y, Htot_RK1_z,
                k1x, k1y, k1z); // update the variables using pointers

            // ---------- RK4 stage 2 ----------
            delm(mx + 0.5*dt*k1x, my + 0.5*dt*k1y, mz + 0.5*dt*k1z,
                Htot_RK2_x, Htot_RK2_y, Htot_RK2_z,
                k2x, k2y, k2z); // update the variables using pointers

            // ---------- RK4 stage 3 ----------
            delm(mx + 0.5*dt*k2x, my + 0.5*dt*k2y, mz + 0.5*dt*k2z,
                Htot_RK3_x, Htot_RK3_y, Htot_RK3_z,
                k3x, k3y, k3z); // update the variables using pointers

            // ---------- RK4 stage 4 ----------
            delm(mx + dt*k3x, my + dt*k3y, mz + dt*k3z,
                Htot_RK4_x, Htot_RK4_y, Htot_RK4_z,
                k4x, k4y, k4z); // update the variables using pointers

            // ---------- Final update ----------
            mx += (dt/6.0)*(k1x + 2.0*k2x + 2.0*k3x + k4x);
            my += (dt/6.0)*(k1y + 2.0*k2y + 2.0*k3y + k4y);
            mz += (dt/6.0)*(k1z + 2.0*k2z + 2.0*k3z + k4z);

            // normalize(mx, my, mz);

            Mx[idx] = mx;
            My[idx] = my;
            Mz[idx] = mz;
            }
        } 
    }
         
    for (double mx : My)    cout << mx << " ";

    ofstream fout("test.dat");
    if (!fout.is_open()) {
        cerr << "Error opening file test.dat!" << endl;
        return 1; // or handle error appropriately
    }
    for (int i = 0; i < 5; ++i) {
        fout << Mx[i] << " " << My[i] << " " << Mz[i] << "\n";
    }
    fout.close();


    return 0;
}
