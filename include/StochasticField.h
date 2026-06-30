#ifndef STOCHASTICFIELD_H
#define STOCHASTICFIELD_H

#include <random>
#include <cmath>

class StochasticField {
public:
    static void generate(
        double D, // noise strength: controls thermal fluctuation amplitude
        double dt,
        std::mt19937& gen,
        double& eta_x,
        double& eta_y,
        double& eta_z
    ) {
        static std::normal_distribution<double> gauss(0.0, 1.0);

        const double sigma = std::sqrt(2.0 * D * dt);

        eta_x = sigma * gauss(gen);
        eta_y = sigma * gauss(gen);
        eta_z = sigma * gauss(gen);
    }
};

#endif