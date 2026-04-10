#include "sofama.h"
#include "agents.h"
#include <cmath>
#include <random>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <vector>

/**
 * 3) F - Целевая функция (Fitness).
 * SOFAMA работает на максимизацию. Если исходная задача — минимизация J (ошибки),
 * мы используем сигмоидальное преобразование, чтобы перевести J в диапазон (0, 1].
 */
double getFitness(double J) {
    return 1.0 / (1.0 + std::exp(J)); // Чем меньше ошибка J, тем ближе Fitness к 0.5-1.0
}

/**
 * Основная функция запуска алгоритма SOFAMA.
 * p - параметры (M, Q, K и т.д.), targetFunc - вычисляемый критерий качества (F).
 */
std::vector<double> runSofama(const SofamaParams& p, std::function<double(const std::vector<double>&)> targetFunc) {
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    
    // Инициализация локальных переменных из структуры параметров
    int N = p.N;                // Размерность (количество коэффициентов 'a')
    double Q = p.Q;              // Порог для выбора между сильной и слабой мутацией
    double gamma_thresh = p.gamma; // Порог отсева по модели Лотки-Вольтерры
    int N_weak = 20000;         // Смещение для слабой мутации (замедляет сходимость)

    // 
    std::vector<std::vector<double>> Z(p.M, std::vector<double>(N)); 
    std::vector<double> J_val(p.M); 
    std::vector<double> F(p.M);     
    
    double best_J = 1e9;            // Переменная для хранения лучшей найденной ошибки
    std::vector<double> best_z(N);  // Вектор лучшего найденного решения

    // 4) Генерация начальной популяции M точек
    for (int i = 0; i < p.M; ++i) {
        for (int j = 0; j < N; ++j) {
            
            std::uniform_real_distribution<double> dist(A_BOUNDS[j].low, A_BOUNDS[j].high);
            Z[i][j] = dist(gen);
        }
        J_val[i] = targetFunc(Z[i]); // Считаем ошибку для созданной точки
        F[i] = getFitness(J_val[i]); // Преобразуем в фитнес для максимизации
        
        // Обновляем глобально лучшее решение, если нашли точку лучше
        if (J_val[i] < best_J) {
            best_J = J_val[i];
            best_z = Z[i];
        }
    }
    
    // 8) k - Счетчик итераций. Начинаем с M+1 (после создания начальной базы)
    int current_pop_size = p.M;
    for (int k = p.M + 1; k <= p.K; ++k) {
        
        // 10) phi(k) - Управляющая последовательность (растет к бесконечности)
        // Регулирует отбор
        double phi_k = std::pow((static_cast<double>(k) / 1000.0 + 1.0), 0.5);
        
        // 11) h(k) 
        double H_k = 0.5 * (1.0 - std::tanh(static_cast<double>(k) - 10000.0));
        
        // zeta_k - Параметр масштаба для распределения Коши (уменьшается со временем)
        double zeta_k = std::pow(1.0 / static_cast<double>(k), 1.0 / (2.0 * N));
        // Для слабой мутации используется большее "время", чтобы шаги были еще меньше
        double zeta_k_weak = std::pow(1.0 / static_cast<double>(k + N_weak), 1.0 / (2.0 * N));

        // 5) Поиск лучшей точки в текущей популяции
        int best_idx = 0;
        double max_F = -1.0;
        std::vector<double> prob(current_pop_size); // Вектор вероятностей выбора агентов
        double sum_prob = 0.0;

        for (int i = 0; i < current_pop_size; ++i) {
            if (F[i] > max_F) {
                max_F = F[i];
                best_idx = i; // Индекс лидера популяции
            }
            // Вероятность выбора родителя пропорциональна его фитнесу в степени phi
            prob[i] = std::pow(F[i], phi_k); 
            sum_prob += prob[i];
        }

        std::vector<double> new_agent(N); 
        std::uniform_real_distribution<double> rand01(0.0, 1.0); 
        
        // Если популяция слишком мала, мутируем просто лучшего агента
        if (current_pop_size < 3) {
            new_agent = Z[best_idx];
            double sigma_weak = zeta_k_weak;
            for (int i = 0; i < N; ++i) {
                double a_L = A_BOUNDS[i].low;
                double a_R = A_BOUNDS[i].high;
                // Генерация новой координаты через обратную функцию распределения Коши
                double arctan_R = std::atan((a_R - new_agent[i]) / sigma_weak);
                double arctan_L = std::atan((a_L - new_agent[i]) / sigma_weak);
                double B_ik = 1.0 / (arctan_R - arctan_L);
                new_agent[i] = new_agent[i] + sigma_weak * std::tan(rand01(gen) / B_ik + arctan_L);
            }
        } else {
            // Выбор 3-х родителей на основе их фитнеса (рулетка)
            std::discrete_distribution<int> select_agent(prob.begin(), prob.end());
            int w_idx = select_agent(gen); // "Центральный" родитель
            int v_idx, u_idx;

            // Гарантируем, что выбрали разных родителей
            int attempts = 0;
            do { v_idx = select_agent(gen); attempts++; } while (v_idx == w_idx && attempts < 100);
            if (v_idx == w_idx) v_idx = (w_idx + 1) % current_pop_size;

            attempts = 0;
            do { u_idx = select_agent(gen); attempts++; } while ((u_idx == v_idx || u_idx == w_idx) && attempts < 100);
            if (u_idx == v_idx || u_idx == w_idx) {
                u_idx = (v_idx + 1) % current_pop_size;
                if (u_idx == w_idx) u_idx = (u_idx + 1) % current_pop_size;
            }
            
            // Анизотропная мутация для каждой координаты i (от a0 до a6)
            for (int i = 0; i < N; ++i) {
                double m_ik, sigma_ik;
                // Сильная/слабая мутация
                if (rand01(gen) < Q) {
                    // Сильная мутация: смещение в сторону лучшего + разность двух случайных
                    m_ik = Z[w_idx][i] + H_k * (Z[best_idx][i] - Z[w_idx][i]) + H_k * (Z[v_idx][i] - Z[u_idx][i]);
                    sigma_ik = zeta_k; // Большой шаг
                } else {
                    // Слабая мутация: малые колебания вокруг выбранного родителя
                    m_ik = Z[w_idx][i];
                    sigma_ik = zeta_k_weak; // Маленький шаг
                }

                // Математический расчет распределения Коши, ограниченного рамками [a_L, a_R]
                double a_L = A_BOUNDS[i].low;
                double a_R = A_BOUNDS[i].high;
                double arctan_R = std::atan((a_R - m_ik) / sigma_ik);
                double arctan_L = std::atan((a_L - m_ik) / sigma_ik);
                double B_ik = 1.0 / (arctan_R - arctan_L);
                
                // Получаем новое значение координаты (ген)
                double mutated_gene = m_ik + sigma_ik * std::tan(rand01(gen) / B_ik + arctan_L);
                
                // 2) П - Жесткий контроль границ (clamping)
                if (mutated_gene > a_R) mutated_gene = a_R;
                if (mutated_gene < a_L) mutated_gene = a_L;
                new_agent[i] = mutated_gene;
            }
        }

        // Оценка нового полученного агента
        double new_J = targetFunc(new_agent); // Считаем его ошибку
        double new_F = getFitness(new_J);     // Считаем его фитнес

        // Если новый агент — лучший за всю историю, сохраняем его
        if (new_J < best_J) {
            best_J = new_J;
            best_z = new_agent;
        }

        // Включение нового агента в текущую популяцию (мю = 1)
        Z.push_back(new_agent);
        F.push_back(new_F);
        current_pop_size++;
        
        // Поиск лидера для гарантированного сохранения при отсеве
        int real_best_idx = 0;
        double current_max_F = -1.0;
        for (int i = 0; i < current_pop_size; ++i) {
            if (F[i] > current_max_F) {
                current_max_F = F[i];
                real_best_idx = i;
            }
        }

        // Нормировочная константа для модели выживания
        double sum_F_phi = 0.0;
        for (double f_val : F) {
            sum_F_phi += std::pow(f_val, phi_k);
        }

        // Временные контейнеры для следующего поколения
        std::vector<std::vector<double>> next_Z;
        std::vector<double> next_F;
        
        // 7) Отсев слабых по модели Лотки-Вольтерры (динамика популяции)
        for (int i = 0; i < current_pop_size; ++i) {
            // Метрика выживания (относительная сила агента)
            double survival_metric = std::pow(F[i], phi_k) / sum_F_phi;
            // Агент выживает, если он сильнее порога гамма ИЛИ если он абсолютный лидер
            if (survival_metric >= gamma_thresh || i == real_best_idx) { 
                next_Z.push_back(Z[i]);
                next_F.push_back(F[i]);
            }
        }
        
        // Обновляем популяцию (удаляем выбывших агентов из памяти)
        Z = std::move(next_Z);
        F = std::move(next_F);
        current_pop_size = Z.size();
    }

    // Возвращаем координаты точки, обеспечившей наименьшую ошибку
    return best_z;
}