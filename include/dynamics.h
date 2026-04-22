#ifndef DYNAMICS_H
#define DYNAMICS_H
#include <vector>

// Задаем константы: аэродинамические коэффициенты БП, целевые параметры и время
struct AeroConstants {
    double k_beta, k_gamma, l_beta, l_gamma, l_psi, l_e, n_beta, n_gamma, n_psi, n_c;
    double alpha_0, omega_y_star, t_T, epsilon;
};

// Результат симуляции - вектор фазовых координат
using SimulationResult = std::vector<double>;

SimulationResult solveSystemODE(const std::vector<double>& a, const AeroConstants& config);


struct TimeStep {
    double t;
    std::vector<double> y;
};
typedef std::vector<TimeStep> Trajectory;

Trajectory getFullTrajectory(const std::vector<double>& a, const AeroConstants& config);

#endif