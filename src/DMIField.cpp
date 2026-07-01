/*
Computes the atomistic DMI contribution to the effective magnetic
field at lattice site (i, j), as a direct lattice sum with a fixed
DMI vector D = (Dx, Dy, Dz).

1. Gathers the four nearest-neighbour spins in the 2D grid and sums
   their components (Sx, Sy, Sz).

2. The Dzyaloshinskii-Moriya interaction couples each spin to its
   nearest neighbors via a cross product, favoring canting / chiral
   textures rather than collinear alignment.

3. Computes the DMI field as the lattice sum
        h_DMI,i = - sum_j D x S_j = - D x (sum of neighbour spins),
   giving (Hdmi_x, Hdmi_y, Hdmi_z) at site (i, j).

4. Combined later with exchange, external, anisotropy, thermal, etc.
   to form the total effective field used in the LLG update.
*/

#include "DMIField.h"
#include "utils.h"

DMIField::DMIField(int Nx, int Ny, double Dx, double Dy, double Dz)
    : Dx(Dx), Dy(Dy), Dz(Dz), Nx_(Nx), Ny_(Ny) {}

// -------------------- Field calculations --------------------
void DMIField::calculate(
    int i, int j,
    double* Mx, double* My, double* Mz,
    double& Hdmi_x, double& Hdmi_y, double& Hdmi_z
    ){
    int ip = (i + 1) % Nx_;        // right neighbor
    int im = (i - 1 + Nx_) % Nx_;  // left  neighbor
    int jp = (j + 1) % Ny_;        // up    neighbor
    int jm = (j - 1 + Ny_) % Ny_;  // down  neighbor

    int ipj  = convert_idx_2dto1D(ip, j,  Nx_, Ny_);
    int imj  = convert_idx_2dto1D(im, j,  Nx_, Ny_);
    int ijp  = convert_idx_2dto1D(i,  jp, Nx_, Ny_);
    int ijm  = convert_idx_2dto1D(i,  jm, Nx_, Ny_);

    // Lattice sum of the four nearest-neighbour spins
    double Sx = Mx[ipj] + Mx[imj] + Mx[ijp] + Mx[ijm];
    double Sy = My[ipj] + My[imj] + My[ijp] + My[ijm];
    double Sz = Mz[ipj] + Mz[imj] + Mz[ijp] + Mz[ijm];

    // h_DMI,i = - ( D x Ssum )
    Hdmi_x = -(Dy * Sz - Dz * Sy);
    Hdmi_y = -(Dz * Sx - Dx * Sz);
    Hdmi_z = -(Dx * Sy - Dy * Sx);
}