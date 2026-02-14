#include "DelmLlg.h"
#include <random>
#include <cmath>




void DelmLlg::calculate(
    // mi's
    double mx, 
    double my, 
    double mz,

    // total external fields
    double Htx, 
    double Hty, 
    double Htz,

    //parameters 
    double randD, 
    double dt, 
    double alpha, 
    double gamma_gyro,

    // RNG
    std::mt19937 &gen,
    std::normal_distribution<double> &gauss,

    // output variables
    double &dmx, double &dmy, double &dmz){

    // random fields
    double eta_x = std::sqrt(2.0 * randD * dt) * gauss(gen);
    double eta_y = std::sqrt(2.0 * randD * dt) * gauss(gen);
    double eta_z = std::sqrt(2.0 * randD * dt) * gauss(gen);

    // Add random field with the total fields from all sources
    Htx += eta_x;
    Hty += eta_y;
    Htz += eta_z;

    // m \times Ht
    double cx = my * Htz - mz * Hty;
    double cy = mz * Htx - mx * Htz;
    double cz = mx * Hty - my * Htx;

    // - m cross (m cross Ht)
    double mdh = mx * Htx + my * Hty + mz * Htz;

    double dcx = Htx - mdh * mx;
    double dcy = Hty - mdh * my;
    double dcz = Htz - mdh * mz;

    double pref = -gamma_gyro / (1.0 + alpha * alpha);

    dmx = pref * (cx + alpha * dcx);
    dmy = pref * (cy + alpha * dcy);
    dmz = pref * (cz + alpha * dcz);
}






