#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <random>
#include <fstream>
#include <string>

#define M_PI 3.14159265358979323846
#include <cmath>

#include "dynamics.h"
#include "criteria.h"
#include "sofama.h"
#include "target_states.h"
#include "agents.h"

// Раскомментируйте, если файл stability.h существует и требуется:
// #include "stability.h" 

using namespace std;

// --- ТЕСТОВЫЕ ФУНКЦИИ ---
double rastrigin(const std::vector<double>& x) {
    double sum = 0;
    for (double val : x) {
        double arg = 2.0 * M_PI * val;
        sum += (val * val - 10.0 * cos(arg) + 10.0);
    }
    return sum;
}

double sphere(const std::vector<double>& x) {
    double sum = 0;
    for (double val : x) sum += val * val;
    return sum;
}

double rosenbrock(const std::vector<double>& x) {
    double sum = 0;
    for (size_t i = 0; i < x.size() - 1; ++i) {
        sum += 100.0 * pow(x[i+1] - x[i]*x[i], 2) + pow(1.0 - x[i], 2);
    }
    return sum;
}

void runConvergenceTests() {
    int N = 5; // Размерность тестового пространства
    // Задаем границы [-5.12, 5.12] классические для Растригина
    std::vector<Bound> test_bounds(N, {-5.12, 5.12}); 
    // struct SofamaParams {
    //     int M;       // M - размер начальной популяции
    //     double Q;    // Q - mutation intensity threshold
    //     int N;       // n - размерность пространства
    //     double gamma;// γ - population membership threshold
    //     int K;       // K – total number of steps
    //     std::vector<Bound> bounds; // Добавлено: границы поиска для универсальности
    // };
    SofamaParams test_p = {500, 0.8, N, 0.005, 500000, test_bounds}; 

    auto rastriginf = [](const std::vector<double>& a) { return rastrigin(a); };
    auto spheref    = [](const std::vector<double>& a) { return sphere(a); };
    auto rosenf     = [](const std::vector<double>& a) { return rosenbrock(a); };

    std::cout << "=== Проверка сходимости SOFAMA ===\n";
    
    std::vector<std::pair<std::string, std::function<double(const std::vector<double>&)>>> tests = {
        {"Rastrigin (min: 0)", rastriginf},
        {"Sphere (min: 0)", spheref},
        {"Rosenbrock (min: 0)", rosenf}
    };

    for(auto& test : tests) {
        auto res = runSofama(test_p, test.second);
        double val = test.second(res);
        std::cout << test.first << " -> Результат: " << std::scientific << val << " (Коорд: ";
        for(double x : res) std::cout << std::fixed << std::setprecision(3) << x << " ";
        std::cout << ")\n";
    }
    std::cout << "==================================\n\n";
}

// --- ЭКСПОРТ ДАННЫХ И АНАЛИЗ УСТОЙЧИВОГО РЕЖИМА ---

// Функция экспорта траекторий в CSV для последующего построения графиков
void exportSimulationData(const std::vector<double>& a, const AeroConstants& c, const std::string& filename) {
    Trajectory traj = getFullTrajectory(a, c);
    
    double i_be = (a[0] - c.l_beta) / c.l_e;
    double i_xe = (a[1] - c.l_gamma) / c.l_e;
    double i_ye = (c.l_psi - a[2]) / c.l_e;
    double q_e  = a[3] / c.l_e;
    double i_bc = (a[4] - c.n_beta) / c.n_c;
    double i_xc = (a[5] - c.n_gamma) / c.n_c;
    double i_yc = (a[6] - c.n_psi) / c.n_c;

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Не удалось открыть файл " << filename << " для записи!\n";
        return;
    }

    out << "t,Y1_beta,Y2_omega_x,Y3_omega_y,Y4_gamma,delta_e,delta_c\n";
    for(auto& step : traj) {
        double delta_e = i_be * step.y[0] + i_xe * step.y[3] + i_ye * step.y[2] + q_e * step.y[1];
        double delta_c = i_bc * step.y[0] + i_xc * step.y[3] + i_yc * step.y[2];
        
        out << step.t << "," << step.y[0] << "," << step.y[1] << "," 
            << step.y[2] << "," << step.y[3] << "," << delta_e << "," << delta_c << "\n";
    }
    out.close();
    std::cout << "Данные фазовых координат экспортированы в " << filename << "\n";
}

// Функция определения характеристик установившегося (устойчивого) режима
void analyzeSteadyState(const std::vector<double>& a, const AeroConstants& config) {
    TargetY ty = calculateTargetStates(a, config);
    Trajectory traj = getFullTrajectory(a, config);
    
    double epsilon = config.epsilon; // константа устойчивого режима (трубка)
    double t_reg = config.t_T; // Время регулирования
    
    // Ищем время вхождения в ε-трубку (с конца траектории) по целевой скорости ω_y
    for (auto it = traj.rbegin(); it != traj.rend(); ++it) {
        if (std::abs(it->y[2] - ty[2]) > epsilon) {
            t_reg = it->t;
            break;
        }
    }
    
    double static_error_omega_y = std::abs(traj.back().y[2] - ty[2]);

    std::cout << "\n--- Характеристики устойчивого режима ---\n";
    std::cout << "Целевая угловая скорость w_y*: " << ty[2] << "\n";
    std::cout << "Фактическая w_y в конце (t=" << config.t_T << "): " << traj.back().y[2] << "\n";
    std::cout << "Статическая ошибка по w_y: " << static_error_omega_y << "\n";
    std::cout << "Время регулирования t_p (при константе ε=" << epsilon << "): " << t_reg << " с\n";
}

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

    // Запуск тестов на глобальный оптимум
    bool test_mode = false;
    if(test_mode) {
        runConvergenceTests();
        return 0;
    } 

    // Данные БПЛА SR200 
    AeroConstants sr200 = {
        0.1946, 0.0883, 47.272, 6.776, 1.742, 176.54, 13.81, 0.108, 0.859, 7.12, 
        3.9, 2.0, 15.0, 0.01 // alpha0, omega_y*, t_T, epsilon (0.01 - константа устойчивости)
    };
    
    // Настройка параметров SOFAMA для задачи синтеза управления
    // Передаем A_BOUNDS из agents.h
    SofamaParams s_p = {50, 0.7, 7, 0.01, 10000, A_BOUNDS}; 
    
    auto func = [&](const std::vector<double>& a) {
        return calculateF(a, sr200);
    };

    std::cout << "Запуск поиска оптимальных коэффициентов управления БПЛА...\n";
    std::vector<double> best_a = runSofama(s_p, func);

    std::cout << "\nОптимальный вектор a: ";
    for(double x : best_a) std::cout << std::fixed << std::setprecision(4) << x << " ";
    
    double final_F = calculateF(best_a, sr200);
    std::cout << "\nИтоговое значение сверточного критерия J: " << final_F << "\n";
    
    printControlGains(best_a, sr200);
    
    // Экспортируем данные для построения графиков
    exportSimulationData(best_a, sr200, "trajectory.csv");


    return 0;
}