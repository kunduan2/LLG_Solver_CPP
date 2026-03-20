/*
===============================================================================
ExchangeField Class
-------------------------------------------------------------------------------
Evaluates the nearest-neighbour exchange field on a two-dimensional
rectangular lattice using a finite-difference Laplacian.

Discrete form (square lattice):

    H_exch(i,j) = C * [
        M(i+1, j) + M(i-1, j)
      + M(i, j+1) + M(i, j-1)
      - 4 M(i, j)
    ]

where:
    C  → exchange prefactor (includes material constants and Δx⁻²)
    4  → number of nearest neighbours on a 2D square lattice

Memory layout:
    Magnetization components are stored in flattened 1D arrays
    using row-major indexing:

        idx(i, j) = i * Nx + j

Boundary treatment:
    Open boundary conditions (no periodic wrapping).
    Spins at the edges interact only with existing neighbours.
    No out-of-bounds memory access occurs.

Implementation goal:
    Provide a stable and well-defined discrete exchange operator
    suitable for LLG time integration.
===============================================================================
*/

#ifndef ExchangeField_H
#define ExchangeField_H

class ExchangeField{
public:
    // Constructor: specify grid size
    ExchangeField(int Nx, int Ny);
    
    // Calculate exchange field at position (i,j)
    void calculate(
        int i, int j,
        double* Mx, double* My, double* Mz,
        double &Hexch_x, double &Hexch_y, double &Hexch_z
    );

private:
    int Nx_, Ny_; // grid dimensions
};


#endif