#ifndef SOFAMA_H
#define SOFAMA_H

#include <vector>
#include <functional>
#include "agents.h"

struct SofamaParams {
    int M;       // M - размер начальной популяции
    double Q;    // Q - mutation intensity threshold
    int N;       // n - размерность пространства
    double gamma;// γ - population membership threshold
    int K;       // K – total number of steps
    std::vector<Bound> bounds; // Добавлено: границы поиска для универсальности
};

double getFitness(double J);
std::vector<double> runSofama(const SofamaParams& p, std::function<double(const std::vector<double>&)> targetFunc);

#endif