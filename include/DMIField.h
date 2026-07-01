#ifndef DMI_FIELD_H
#define DMI_FIELD_H

/*
DMIField  (atomistic / lattice form, fixed D vector)

Computes the Dzyaloshinskii-Moriya contribution to the effective field
at lattice site (i, j) on a 2D grid, as a direct lattice sum over the
four nearest neighbours using a single fixed DMI vector D = (Dx,Dy,Dz).

Energy:  F_DMI = sum_<ij> D . (S_i x S_j)
Field :  h_DMI,i = - sum_j  D x S_j
                 = - D x ( S_ip + S_im + S_ijp + S_ijm )

with j running over the four nearest neighbours. D is the same fixed
vector on every bond and is supplied as input (Dx, Dy, Dz). The lattice
constant is absorbed into the components of D.

Component form of  - ( D x Ssum ),  Ssum = (Sx, Sy, Sz):
    h_x = - ( Dy*Sz - Dz*Sy )
    h_y = - ( Dz*Sx - Dx*Sz )
    h_z = - ( Dx*Sy - Dy*Sx )
*/

class DMIField {
public:
    DMIField(int Nx, int Ny, double Dx, double Dy, double Dz);

    void calculate(
        int i, int j,
        double* Mx, double* My, double* Mz,
        double& Hdmi_x, double& Hdmi_y, double& Hdmi_z
    );

    double Dx;   // DMI vector x-component (fixed, lattice spacing absorbed)
    double Dy;   // DMI vector y-component
    double Dz;   // DMI vector z-component

private:
    int Nx_;
    int Ny_;
};

#endif // DMI_FIELD_H