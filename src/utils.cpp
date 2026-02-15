// #include "utils.h"
// #include <cmath>
// #include <random>

// // -------------------- RNG (translation-unit local) --------------------
// static std::mt19937 gen(42);
// static std::normal_distribution<double> gauss(0.0, 1.0);

// // -------------------- Utility functions --------------------
// void normalize(double &mx, double &my, double &mz)
// {
//     double n = std::sqrt(mx*mx + my*my + mz*mz);
//     mx /= n;
//     my /= n;
//     mz /= n;
// }

// int convert_idx_2dto1D(int i, int j, int Nx, int Ny)
// {
//     return i * Nx + j;
// }

// int indices_conversion_3dto1D(int i, int j, int k,
//                               int Nx, int Ny, int Nz)
// {
//     return k * Nx * Ny + i * Nx + j;
// }

// int pbc(int i, int N)
// {
//     return (i + N) % N;
// }

// double mean(double* A, int size)
// {
//     double s = 0.0;
//     for (int i = 0; i < size; i++)
//         s += A[i];
//     return s / size;
// }

// -------------------- Field calculations --------------------
// void exchange_field(int i, int j, int Nx, int Ny,
//                     double* Mx, double* My, double* Mz,
//                     double &Hexch_x, double &Hexch_y, double &Hexch_z)
// {
//     int ip = (i + 1) % Nx;
//     int im = (i - 1 + Nx) % Nx;
//     int jp = (j + 1) % Ny;
//     int jm = (j - 1 + Ny) % Ny;

//     int ij   = convert_idx_2dto1D(i,  j,  Nx, Ny);
//     int ipj  = convert_idx_2dto1D(ip, j,  Nx, Ny);
//     int imj  = convert_idx_2dto1D(im, j,  Nx, Ny);
//     int ijp  = convert_idx_2dto1D(i,  jp, Nx, Ny);
//     int ijm  = convert_idx_2dto1D(i,  jm, Nx, Ny);

//     Hexch_x = Mx[ipj] + Mx[imj] + Mx[ijp] + Mx[ijm] - 4.0 * Mx[ij];
//     Hexch_y = My[ipj] + My[imj] + My[ijp] + My[ijm] - 4.0 * My[ij];
//     Hexch_z = Mz[ipj] + Mz[imj] + Mz[ijp] + Mz[ijm] - 4.0 * Mz[ij];
// }

// void H_of_t(double t, double t_switch, double Bext,
//             double &Hx, double &Hy, double &Hz){
//     Hx = 0.0;
//     Hy = 0.0;
//     Hz = -Bext;
// }

// // -------------------- LLG derivative --------------------
// void delm(double mx, double my, double mz,
//           double Hx, double Hy, double Hz,
//           double D, double dt, double alpha, double gamma_gyro,
//           std::mt19937 &gen,
//           std::normal_distribution<double> &gauss,
//           double &dmx, double &dmy, double &dmz){
                        
//     double eta_x = std::sqrt(2.0 * D * dt) * gauss(gen);
//     double eta_y = std::sqrt(2.0 * D * dt) * gauss(gen);
//     double eta_z = std::sqrt(2.0 * D * dt) * gauss(gen);

//     double Htx = Hx + eta_x;
//     double Hty = Hy + eta_y;
//     double Htz = Hz + eta_z;

//     double cx = my * Htz - mz * Hty;
//     double cy = mz * Htx - mx * Htz;
//     double cz = mx * Hty - my * Htx;

//     double mdh = mx * Htx + my * Hty + mz * Htz;

//     double dcx = Htx - mdh * mx;
//     double dcy = Hty - mdh * my;
//     double dcz = Htz - mdh * mz;

//     double pref = -gamma_gyro / (1.0 + alpha * alpha);

//     dmx = pref * (cx + alpha * dcx);
//     dmy = pref * (cy + alpha * dcy);
//     dmz = pref * (cz + alpha * dcz);
// }

// -------------------- Solver --------------------
// void llg_solver(double t_switch, double Bext, double Aexch,
//                 int Nt, int Nx, int Ny,
//                 double t, double dt, double D,
//                 double alpha, double gamma_gyro,
//                 double* Mx, double* My, double* Mz)
// {
//     double Hex, Hey, Hez;
//     double H1x, H1y, H1z, H2x, H2y, H2z, H3x, H3y, H3z, H4x, H4y, H4z;
//     double k1x, k1y, k1z, k2x, k2y, k2z, k3x, k3y, k3z, k4x, k4y, k4z;

//     for (int tn = 1; tn < Nt; tn++) {
//         for (int i = 0; i < Ny - 1; i++) {
//             for (int j = 0; j < Nx - 1; j++) {

//                 exchange_field(i, j, Nx, Ny, Mx, My, Mz, Hex, Hey, Hez);

//                 H_of_t(t,           t_switch, Bext, H1x, H1y, H1z);
//                 H_of_t(t + 0.5*dt,  t_switch, Bext, H2x, H2y, H2z);
//                 H_of_t(t + 0.5*dt,  t_switch, Bext, H3x, H3y, H3z);
//                 H_of_t(t + dt,      t_switch, Bext, H4x, H4y, H4z);

//                 int idx = convert_idx_2dto1D(i, j, Nx, Ny);

//                 double mx = Mx[idx];
//                 double my = My[idx];
//                 double mz = Mz[idx];

//                 delm(mx, my, mz, H1x + Aexch*Hex, H1y + Aexch*Hey, H1z + Aexch*Hez,
//                      D, dt, alpha, gamma_gyro, gen, gauss, k1x, k1y, k1z);

//                 delm(mx + 0.5*dt*k1x, my + 0.5*dt*k1y, mz + 0.5*dt*k1z,
//                      H2x + Aexch*Hex, H2y + Aexch*Hey, H2z + Aexch*Hez,
//                      D, dt, alpha, gamma_gyro, gen, gauss, k2x, k2y, k2z);

//                 delm(mx + 0.5*dt*k2x, my + 0.5*dt*k2y, mz + 0.5*dt*k2z,
//                      H3x + Aexch*Hex, H3y + Aexch*Hey, H3z + Aexch*Hez,
//                      D, dt, alpha, gamma_gyro, gen, gauss, k3x, k3y, k3z);

//                 delm(mx + dt*k3x, my + dt*k3y, mz + dt*k3z,
//                      H4x + Aexch*Hex, H4y + Aexch*Hey, H4z + Aexch*Hez,
//                      D, dt, alpha, gamma_gyro, gen, gauss, k4x, k4y, k4z);

//                 mx += (dt/6.0)*(k1x + 2*k2x + 2*k3x + k4x);
//                 my += (dt/6.0)*(k1y + 2*k2y + 2*k3y + k4y);
//                 mz += (dt/6.0)*(k1z + 2*k2z + 2*k3z + k4z);

//                 normalize(mx, my, mz);

//                 Mx[idx] = mx;
//                 My[idx] = my;
//                 Mz[idx] = mz;
//             }
//         }
//         t += dt;
//     }
// }
