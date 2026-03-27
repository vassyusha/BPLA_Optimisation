#include "target_states.h"
#include <vector>

// Вычисление устойчивых состояний Y* согласно формулам (27-31)
TargetY calculateTargetStates(const std::vector<double>& a, const AeroConstants& config) {
    // Вектор для хранения Y1*, Y2*, Y3*, Y4*, Y6*, (Y8* = 0)
    TargetY ty(6); 

    double w_y_star = config.omega_y_star;
    double k_beta = config.k_beta;
    double k_gamma = config.k_gamma;
    double l_psi = config.l_psi;

    // Y1* = - (a6 * w_y*) / a4  (27)
    ty[0] = -(a[6] * w_y_star) / a[4];

    // Y2* = 0  (28)
    ty[1] = 0.0;

    // Y3* = w_y* (29)
    ty[2] = w_y_star;

    // Y4* = - (k_beta * a6 * w_y*) / a4*k_gamma - l_psi*w_y*/a3  (30)
    ty[3] = -((k_beta * a[6] * w_y_star)/a[4] + w_y_star) / k_gamma; 

    // Y6* = (a0 * a6 * w_y*) / a4 + (l_psi * w_y*) / a3   (31)
    ty[4] = ( (a[0] * a[6] * w_y_star / a[4]) + l_psi*w_y_star)/a[3];

    // Y8* = 0
    ty[5] = 0;
    return ty;
}