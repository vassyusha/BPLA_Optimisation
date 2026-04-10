#include "dynamics.h"
#include <vector>
#include <cmath>

using namespace std;

// Решаем систему уравнений (21)-(25) из постановки
vector<double> derivatives(const vector<double>& y, const vector<double>& a, const AeroConstants& config) {
    vector<double> dy(6);

    // Уравнение (21): Y1 = -k_beta*Y1 + alpha_0*Y2 + Y3 + k_gamma*Y4
    dy[0] = -config.k_beta * y[0] + config.alpha_0 * y[1] + y[2] + config.k_gamma * y[3];

    // Уравнение (22): Y2 = -a0*Y1 - a1*Y2 - a2*Y3 - w_y* + a3*Y6 - l_psi*w_y*
    dy[1] = -a[0] * y[0] - a[1] * y[1] - a[2] * y[2] - config.omega_y_star + a[3] * y[4] - config.l_psi * config.omega_y_star;

    // Уравнение (23): Y3 = -a4*Y1 - a5*Y2 - a6*Y3
    dy[2] = -a[4] * y[0] - a[5] * y[1] - a[6] * y[2];

    // Уравнение (24): Y4 = Y2
    dy[3] = y[1];

    // Уравнение (25): Y6 = Y3 - w_y*
    dy[4] = y[2] - config.omega_y_star;

    // Уравнение для интеграла (Y3 - w_y*)^2 (свертка)
    double diff = y[2] - config.omega_y_star;
    dy[5] = diff * diff;

    return dy;
}

SimulationResult solveSystemODE(const vector<double>& a, const AeroConstants& config) {
    // Начальные условия (состояние покоя в начале маневра)
    vector<double> y = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    double t = 0.0;
    double dt = 0.01; // Шаг интегрирования
    double t_T = config.t_T; // Конечное время разворота

    // Метод Рунге-Кутты 4-го порядка
    while (t < t_T) {
        if (t + dt > t_T) dt = t_T - t; 

        vector<double> k1 = derivatives(y, a, config);
        
        vector<double> y_k2(6);
        for(int i=0; i<6; ++i) y_k2[i] = y[i] + 0.5 * dt * k1[i];
        vector<double> k2 = derivatives(y_k2, a, config);

        vector<double> y_k3(6);
        for(int i=0; i<6; ++i) y_k3[i] = y[i] + 0.5 * dt * k2[i];
        vector<double> k3 = derivatives(y_k3, a, config);

        vector<double> y_k4(6);
        for(int i=0; i<6; ++i) y_k4[i] = y[i] + dt * k3[i];
        vector<double> k4 = derivatives(y_k4, a, config);

        for (int i = 0; i < 6; ++i) {
            y[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        t += dt;
    }
    return y; 
}