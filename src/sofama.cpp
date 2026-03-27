#include "sofama.h"
#include "agents.h"
#include <cmath>
#include <random>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <vector>

// Трансформация целевой функции: минимизация J -> максимизация F 
double getFitness(double J) {
    return 1.0 / (1.0 + std::exp(J)); // J = наш сверточный критерий
}

std::vector<double> runSofama(const SofamaParams& p, std::function<double(const std::vector<double>&)> targetFunc) {
    std::random_device rd;
    std::mt19937 gen(rd());
    int N = 7; // Размерность пространства (a0...a6)
    
    // Параметры алгоритма SoFAMA (можно вынести в SofamaParams)
    double Q = 0.7; // Порог интенсивности мутации (G в некоторых частях текста) [cite: 988]
    int N_weak = 20000; // Сдвиг шага для слабой мутации [cite: 988]
    double gamma_thresh = 0.01; // Порог принадлежности к популяции [cite: 988]

    // 1. Инициализация популяции [cite: 679]
    std::vector<std::vector<double>> Z(p.M, std::vector<double>(N));
    std::vector<double> J_val(p.M);
    std::vector<double> F(p.M);

    double best_J = 1e9;
    std::vector<double> best_z(N);

    for (int i = 0; i < p.M; ++i) {
        for (int j = 0; j < N; ++j) {
            std::uniform_real_distribution<double> dist(A_BOUNDS[j].low, A_BOUNDS[j].high);
            Z[i][j] = dist(gen);
        }
        J_val[i] = targetFunc(Z[i]);
        F[i] = getFitness(J_val[i]);
        
        if (J_val[i] < best_J) {
            best_J = J_val[i];
            best_z = Z[i];
        }
    }

    // Главный цикл алгоритма (начинается с k = m + 1) [cite: 683]
    int current_pop_size = p.M;
    for (int k = p.M + 1; k <= p.K; ++k) {
        if (k % 1000 == 0) {
            std::cout << "Iteration: " << k << " | Pop Size: " << current_pop_size 
              << " | Best J: " << best_J << std::endl;
        }

        // Управляющие последовательности [cite: 850, 988]
        double phi_k = std::pow((static_cast<double>(k) / 1000.0 + 1.0), 0.5);
        double H_k = 0.5 * (1.0 - std::tanh(static_cast<double>(k) - 10000.0));
        double zeta_k = std::pow(1.0 / static_cast<double>(k), 1.0 / (2.0 * N));
        double zeta_k_weak = std::pow(1.0 / static_cast<double>(k + N_weak), 1.0 / (2.0 * N));

        // Выбор лучшего агента текущей популяции [cite: 684]
        int best_idx = 0;
        double max_F = -1.0;
        std::vector<double> prob(current_pop_size);
        double sum_prob = 0.0;

        for (int i = 0; i < current_pop_size; ++i) {
            if (F[i] > max_F) {
                max_F = F[i];
                best_idx = i;
            }
            prob[i] = std::pow(F[i], phi_k); // [cite: 689]
            sum_prob += prob[i];
        }

        // [ПРАВКА 2] Выбор опорных векторов (Статья, раздел "Anisotropic Mutation")
        // Необходимо гарантировать w != v != u для предотвращения вырождения мутации

        // 1. Объявляем переменную ЗАРАНЕЕ
        std::vector<double> new_agent(N);
        std::uniform_real_distribution<double> rand01(0.0, 1.0);
        // 2. Проверяем размер популяции ПЕРЕД выбором индексов
        if (current_pop_size < 3) {
            // Если выродились — просто берем лучшего и чуть-чуть «шатаем» его (слабая мутация)
            new_agent = Z[best_idx];
            double sigma_weak = zeta_k_weak;
            for (int i = 0; i < N; ++i) {
                double a_L = A_BOUNDS[i].low;
                double a_R = A_BOUNDS[i].high;
                double arctan_R = std::atan((a_R - new_agent[i]) / sigma_weak);
                double arctan_L = std::atan((a_L - new_agent[i]) / sigma_weak);
                double B_ik = 1.0 / (arctan_R - arctan_L);
                new_agent[i] = new_agent[i] + sigma_weak * std::tan(rand01(gen) / B_ik + arctan_L);
            }
        } else {
            // Если агентов >= 3, выполняем стандартную SOFAMA мутацию
            std::discrete_distribution<int> select_agent(prob.begin(), prob.end());
            
            int w_idx = select_agent(gen);
            int v_idx, u_idx;
            do { v_idx = select_agent(gen); } while (v_idx == w_idx);
            do { u_idx = select_agent(gen); } while (u_idx == v_idx || u_idx == w_idx);

            for (int i = 0; i < N; ++i) {
                double m_ik, sigma_ik;
                if (rand01(gen) < Q) {
                    m_ik = Z[w_idx][i] + H_k * (Z[best_idx][i] - Z[w_idx][i]) + H_k * (Z[v_idx][i] - Z[u_idx][i]);
                    sigma_ik = zeta_k;
                } else {
                    m_ik = Z[w_idx][i];
                    sigma_ik = zeta_k_weak;
                }

                // Ограниченная мутация Коши (твой код...)
                double a_L = A_BOUNDS[i].low;
                double a_R = A_BOUNDS[i].high;
                double arctan_R = std::atan((a_R - m_ik) / sigma_ik);
                double arctan_L = std::atan((a_L - m_ik) / sigma_ik);
                double B_ik = 1.0 / (arctan_R - arctan_L);
                
                double mutated_gene = m_ik + sigma_ik * std::tan(rand01(gen) / B_ik + arctan_L);
                
                if (mutated_gene > a_R) mutated_gene = a_R;
                if (mutated_gene < a_L) mutated_gene = a_L;
                new_agent[i] = mutated_gene;
            }
        }

        // Вычисляем значение целевой функции для нового агента [cite: 713]
        double new_J = targetFunc(new_agent);
        double new_F = getFitness(new_J);

        // Обновляем лучший глобальный результат
        if (new_J < best_J) {
            best_J = new_J;
            best_z = new_agent;
        }

        // 6. Динамическое обновление популяции по модели Лотки-Вольтерры [cite: 714, 715]
        Z.push_back(new_agent);
        F.push_back(new_F);
        current_pop_size++;
        
        int real_best_idx = 0;
        double current_max_F = -1.0;
        for (int i = 0; i < current_pop_size; ++i) {
            if (F[i] > current_max_F) {
                current_max_F = F[i];
                real_best_idx = i;
            }
        }

        // Вычисляем сумму для нормировки выживания
        double sum_F_phi = 0.0;
        for (double f_val : F) {
            sum_F_phi += std::pow(f_val, phi_k);
        }

        std::vector<std::vector<double>> next_Z;
        std::vector<double> next_F;
        
        for (int i = 0; i < current_pop_size; ++i) {
            double survival_metric = std::pow(F[i], phi_k) / sum_F_phi;
            // Используем обновленный real_best_idx
            if (survival_metric >= gamma_thresh || i == real_best_idx) { 
                next_Z.push_back(Z[i]);
                next_F.push_back(F[i]);
            }
        }
        
        Z = std::move(next_Z);
        F = std::move(next_F);
        current_pop_size = Z.size();
    }

    return best_z;
}