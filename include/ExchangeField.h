#ifndef ExchangeField_H
#define ExchangeField_H

class ExchangeField{
    public:
        // Constructor: specify grid size
        ExchangeField(int Nx, int Ny);
        
        // Calculate exchange field at position (i,j)
        void calculate(int i, int j,
                       double* Mx, double* My, double* Mz,
                       double &Hexch_x, double &Hexch_y, double &Hexch_z);
    private:
        int Nx_, Ny_; // grid dimensions

};


#endif