#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <random>

#include "dynamics.h"
#include "criteria.h"
#include "sofama.h"
#include "target_states.h"
#include "agents.h"
#include "stability.h"

using namespace std;



// функции для проверки работоспособности отдельных функций (это временно, как только всё заработает - можно убрать)
void generate_and_check_stability(); // проверка работоспособности stability.cpp


// Функция перевода в коэффициенты рулей (зависимость в конце файла)
void printControlGains(const std::vector<double>& a, const AeroConstants& c) {
    double i_beta_e = (a[0] - c.l_beta) / c.l_e;
    double i_x_e    = (a[1] - c.l_gamma) / c.l_e;
    double i_y_e    = (c.l_psi - a[2]) / c.l_e;
    double q_e      = a[3] / c.l_e;
    double i_beta_c = (a[4] - c.n_beta) / c.n_c;
    double i_x_c    = (a[5] - c.n_gamma) / c.n_c;
    double i_y_c    = (a[6] - c.n_psi) / c.n_c;

    std::cout << "\n--- Передаточные числа автопилота ---\n";
    std::cout << "Элероны: i_be=" << i_beta_e << ", i_xe=" << i_x_e << ", i_ye=" << i_y_e << ", qe=" << q_e << "\n";
    std::cout << "РН:      i_bc=" << i_beta_c << ", i_xc=" << i_x_c << ", i_yc=" << i_y_c << "\n";
}

int main() {
    system("chcp 65001 > nul");

    // Данные SR200 из файла
    AeroConstants sr200 = {
        0.1946, 0.0883, 47.272, 6.776, 1.742, 176.54, 13.81, 0.108, 0.859, 7.12, // Аэродинамика
        3.9, 2.0, 15.0, 0.5 // alpha0, omega_y*, t_T, epsilon
    };

    // 2. Настройка параметров SOFAMA
    // M = 50 (агентов)
    // Q = 0.3 (тяга к лучшему решению)
    // N = 7 (размерность вектора a: a0...a6)
    // gamma = 1.0 (масштаб Коши; 0.01 был слишком мал для начала поиска)
    // K = 2000 (итераций)
    SofamaParams s_p = {50, 0.3, 7, 1.0, 2000};
    
    auto func = [&](const std::vector<double>& a) {
        return calculateF(a, sr200);
    };

    std::vector<double> best_a = runSofama(s_p, func);

    std::cout << "Оптимальный вектор a: ";
    for(double x : best_a) std::cout << std::fixed << std::setprecision(4) << x << " ";
    
    std::cout << "\nЗначение критерия F: " << calculateF(best_a, sr200);
    
    printControlGains(best_a, sr200);

    return 0;
}





















// тест работоспособности проверки устойчивости
void generate_and_check_stability() {
    // 1. SR200
    AeroConstants sr200 = {
        0.1946, 0.0883, 47.272, 6.776, 1.742, 176.54, 13.81, 0.108, 0.859, 7.12, // Аэродинамика
        3.9, 2.0, 15.0, 0.5 // alpha0, omega_y*, t_T, epsilon
    };

    // 2. Настройка генератора случайных чисел
    mt19937 gen(static_cast<unsigned>(chrono::system_clock::now().time_since_epoch().count()));
    
    // Подготавливаем распределения для каждого a_i согласно вашим интервалам
    // Используем A_BOUNDS из agents.h для соответствия
    vector<uniform_real_distribution<double>> dists;
    for (const auto& b : A_BOUNDS) {
        dists.emplace_back(b.low, b.high);
    }

    cout << "--- Поиск случайных устойчивых точек ---\n";
    int found = 0;
    int attempts = 0;
    while (found < 5 && attempts < 100000) {
        attempts++;
        vector<double> random_a(7);
        for(int i = 0; i < 7; ++i) random_a[i] = dists[i](gen);

        if (isSystemStable(random_a, sr200)) {
            found++;
            cout << "[Устойчив] Точка #" << found << " (попытка " << attempts << ")" << endl;
        }
    }

}