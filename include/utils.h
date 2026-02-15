#ifndef LLG_UTILS_H
#define LLG_UTILS_H

#include <random>
#include <cmath>
#include <iostream>


/*
Functions for general mathematical operations and data manipulation
*/

// RNG (translation-unit local)
static std::mt19937 gen(42);
static std::normal_distribution<double> gauss(0.0, 1.0);

// Normalizes a 3D vector to unit length
inline void normalize(double &mx, double &my, double &mz){
    double n = std::sqrt(mx*mx + my*my + mz*mz);
    mx /= n;
    my /= n;
    mz /= n;
}

inline double norm(double x, double y, double z){
    return x*x + y*y + z*z;

}
// Converts 2D grid indices (i,j) to 1D array index (row-major order)
inline int convert_idx_2dto1D(int i, int j, int Nx, int Ny){
    return i * Nx + j;
}

// Converts 3D grid indices (i,j,k) to 1D array index (row-major order)
inline int indices_conversion_3dto1D(
    int i, int j, int k,
    int Nx, int Ny, int Nz
){
    return k * Nx * Ny + i * Nx + j;
}

// Applies periodic boundary conditions: \equiv: int i_right  = (i == N-1) ? 0 : i+1;
inline int pbc(int i, int N){
    return (i + N) % N;
}

// Computes arithmetic mean of an array of given size
inline double mean(double* A, int size)
{
    double s = 0.0;
    for (int i = 0; i < size; i++)
        s += A[i];
    return s / size;
}

inline void print2Darr(double* A, int rows, int cols){  // n is size
    for(int i=0; i<rows; i++ ){
        for(int j=0; j<cols; j++){
            std::cout<<A[i*rows + j] << " ";
        }
        std::cout<<std::endl;
    }
}

#endif // LLG_UTILS_H
