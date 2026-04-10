#include "target_states.h"
#include <vector>

// Вычисляем целевые значения фазовых координат в конечный момент времени (Формулы 27-31)
TargetY calculateTargetStates(const std::vector<double>& a, const AeroConstants& config) {
    TargetY ty(6); 

    double w_y_star = config.omega_y_star;
    double k_beta = config.k_beta;
    double k_gamma = config.k_gamma;
    double l_psi = config.l_psi;

    // ПРИМЕЧАНИЕ: В тексте статьи формулы выглядят как "Y1*=-a6ωy*a4", без знака дроби.
    // Однако с точки зрения физики аэродинамики деление здесь уместнее. Оставлено деление.
    
    // (27)
    ty[0] = -(a[6] * w_y_star) / a[4];
    // (28)
    ty[1] = 0.0;
    // (29)
    ty[2] = w_y_star;
    // (30)
    ty[3] = -((k_beta * a[6] * w_y_star)/a[4] + w_y_star) / k_gamma; 
    // (31)
    ty[4] = ( (a[0] * a[6] * w_y_star / a[4]) + l_psi*w_y_star)/a[3];
    // Интеграл = 0
    ty[5] = 0;
    
    return ty;
}