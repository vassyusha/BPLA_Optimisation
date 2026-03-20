
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <chrono>
#include <random>
#include <omp.h>

using namespace std;

struct Params {
    double k_beta = 0.1946;
    double k_gamma = 0.0883;
    double alpha_0 = 3.9;
    double l_psi = 1.742;
    double omega_y_star = 2.0; // Согласно алгоритму "задача ЛА"
    double t_T = 15.0;

    // Веса xi согласно пропорциям 1/3.5 : 1/15 : 1/3 : 1/25 : 1/0.05 : 1/0.05
    vector<double> xi;
    double eps = 0.5;

    Params() {
        vector<double> p = { 1.0 / 3.5, 1.0 / 15.0, 1.0 / 3.0, 1.0 / 25.0, 1.0 / 0.05, 1.0 / 0.05 };
        double sum_p = accumulate(p.begin(), p.end(), 0.0);
        for (double v : p) xi.push_back(v / sum_p);
    }
};

double objective_function(const vector<double>& a, const Params& p) {
    // 1. Расчет устойчивого состояния Y* (вывод из системы 21-25)
    if (abs(a[4]) < 1e-6 || abs(a[3]) < 1e-6) return 1e12;

    double Y1_star = -(a[6] * p.omega_y_star) / a[4]; // Из dY3/dt=0
    double Y2_star = 0.0;
    double Y3_star = p.omega_y_star;
    double Y4_star = (p.k_beta * Y1_star - p.omega_y_star) / p.k_gamma; // Из dY1/dt=0
    double Y6_star = (a[0] * Y1_star + p.l_psi * p.omega_y_star) / a[3]; // Из dY2/dt=0

    // 2. Решение системы ОДУ (y0 = {0,0,0,0,0,0})
    vector<double> y = { 0, 0, 0, 0, 0, 0 };
    double dt = 0.01;

    for (double t = 0; t < p.t_T; t += dt) {
        double dy[6];
        // Y1: beta, Y2: omega_x, Y3: omega_y, Y4: gamma, Y5: rho (integral), Y6: Y8 (performance)
        dy[0] = -p.k_beta * y[0] + p.alpha_0 * y[1] + y[2] + p.k_gamma * y[3];
        dy[1] = -a[0] * y[0] - a[1] * y[1] - a[2] * (y[2] - p.omega_y_star) + a[3] * y[4] - p.l_psi * p.omega_y_star;
        dy[2] = -a[4] * y[0] - a[5] * y[1] - a[6] * y[2];
        dy[3] = y[1];
        dy[4] = y[2] - p.omega_y_star;
        dy[5] = pow(y[2] - p.omega_y_star, 2);

        for (int i = 0; i < 6; ++i) y[i] += dy[i] * dt;

        if (std::isnan(y[0]) || abs(y[0]) > 50.0) return 1e12;
    }

    // 3. Частные критерии F1-F6
    double F[6];
    F[0] = abs(y[0] - Y1_star);
    F[1] = abs(y[1]);
    F[2] = abs(y[2] - p.omega_y_star);
    F[3] = abs(y[3] - Y4_star);
    F[4] = abs(y[4] - Y6_star);
    F[5] = y[5]; // Y8(15)

    // 4. Ограничение и свертка
    if (F[5] > p.eps) return 1e10 + F[5];

    double sum_F = 0;
    for (int i = 0; i < 6; ++i) sum_F += p.xi[i] * F[i];
    return sum_F;
}

struct Trial { vector<double> a; double z; };

void run_sofama() {
    auto start_time = chrono::steady_clock::now();
    Params p;
    const int dim = 7;
    // Границы оптимизации для вектора a из "задача ЛА"
    vector<pair<double, double>> bounds = {
        {0, 2}, {4, 12}, {-200, 1.74}, {0, 176.54}, {8, 20}, {0, 3}, {0, 8}
    };

    vector<Trial> trials;
    cout << "SOFAMA: Phase 1 - Global Exploration..." << endl;
#pragma omp parallel
    {
        mt19937 local_gen(omp_get_thread_num() + (int)time(0));
        vector<Trial> local_batch;
#pragma omp for
        for (int i = 0; i < 15000; ++i) {
            vector<double> a(dim);
            for (int j = 0; j < dim; ++j)
                a[j] = uniform_real_distribution<double>(bounds[j].first, bounds[j].second)(local_gen);
            local_batch.push_back({ a, objective_function(a, p) });
        }
#pragma omp critical
        trials.insert(trials.end(), local_batch.begin(), local_batch.end());
    }

    sort(trials.begin(), trials.end(), [](const Trial& a, const Trial& b) { return a.z < b.z; });

    cout << "SOFAMA: Phase 2 - Iterative Refinement..." << endl;
    int iter = 0;
    while (chrono::duration<double>(chrono::steady_clock::now() - start_time).count() < 14.5) {
        vector<Trial> next_generation;
#pragma omp parallel
        {
            mt19937 local_gen(iter + omp_get_thread_num() + (int)time(0));
            vector<Trial> local_results;
#pragma omp for
            for (int i = 0; i < 1000; ++i) {
                int idx = uniform_int_distribution<int>(0, 200)(local_gen);
                vector<double> a = trials[idx].a;
                double radius = 0.1 * exp(-(double)iter / 100.0);
                for (int j = 0; j < dim; ++j) {
                    a[j] += normal_distribution<double>(0, radius * (bounds[j].second - bounds[j].first))(local_gen);
                    a[j] = max(bounds[j].first, min(bounds[j].second, a[j]));
                }
                local_results.push_back({ a, objective_function(a, p) });
            }
#pragma omp critical
            next_generation.insert(next_generation.end(), local_results.begin(), local_results.end());
        }
        trials.insert(trials.end(), next_generation.begin(), next_generation.end());
        sort(trials.begin(), trials.end(), [](const Trial& a, const Trial& b) { return a.z < b.z; });
        if (trials.size() > 5000) trials.resize(1000);
        if (iter % 20 == 0) cout << "Iteration " << iter << " | Best F: " << trials[0].z << endl;
        iter++;
    }

    cout << "\nOPTIMIZATION COMPLETE" << endl;
    cout << "Final F: " << trials[0].z << endl;
    for (int i = 0; i < dim; ++i) printf("a%d = %10.6f\n", i, trials[0].a[i]);
}

int main() {
    run_sofama();
    return 0;
}
