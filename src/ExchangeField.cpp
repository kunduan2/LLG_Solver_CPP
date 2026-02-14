#include "ExchangeField.h"
#include "utils.h"
#include <cmath>

ExchangeField::ExchangeField(int Nx, int Ny): Nx_(Nx), Ny_(Ny) {}

// -------------------- Field calculations --------------------
void ExchangeField::calculate(int i, int j,
                    double* Mx, double* My, double* Mz,
                    double &Hexch_x, double &Hexch_y, double &Hexch_z)
{
    int ip = (i + 1) % Nx_;       // \equiv: int i_right  = (i == N-1) ? 0 : i+1;
    int im = (i - 1 + Nx_) % Nx_;  // \equiv: int i_left  = (i == 0) ? N-1 : i-1;
    int jp = (j + 1) % Ny_;       
    int jm = (j - 1 + Ny_) % Ny_;

    int ij   = convert_idx_2dto1D(i,  j,  Nx_, Ny_);
    int ipj  = convert_idx_2dto1D(ip, j,  Nx_, Ny_);
    int imj  = convert_idx_2dto1D(im, j,  Nx_, Ny_);
    int ijp  = convert_idx_2dto1D(i,  jp, Nx_, Ny_);
    int ijm  = convert_idx_2dto1D(i,  jm, Nx_, Ny_);

    Hexch_x = Mx[ipj] + Mx[imj] + Mx[ijp] + Mx[ijm] - 4.0 * Mx[ij];
    Hexch_y = My[ipj] + My[imj] + My[ijp] + My[ijm] - 4.0 * My[ij];
    Hexch_z = Mz[ipj] + Mz[imj] + Mz[ijp] + Mz[ijm] - 4.0 * Mz[ij];
}





/*
Note:
0 % 5 = 0
1 % 5 = 1
2 % 5 = 2
3 % 5 = 3
4 % 5 = 4
5 % 5 = 0
6 % 5 = 1

*/