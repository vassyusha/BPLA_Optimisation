#ifndef SOFAMA_H
#define SOFAMA_H
#include <vector>
#include <functional>

struct SofamaParams {
    int M;       // 4) M - размер начальной популяции
    double Q;    // 5) Q (или G) – mutation intensity threshold
    int N;       // 1) n - размерность пространства
    double gamma;// 7) γ - population membership threshold
    int K;       // 9) K – total number of steps
};


double getFitness(double J);

std::vector<double> runSofama(const SofamaParams& p, std::function<double(const std::vector<double>&)> targetFunc);
#endif