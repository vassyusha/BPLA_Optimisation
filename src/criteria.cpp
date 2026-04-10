#include "criteria.h"
#include "target_states.h"
#include "dynamics.h"
#include <cmath>

double calculateF(const std::vector<double>& a, const AeroConstants& config) {
    // Получаем текущие значения в момент t_T
    SimulationResult y = solveSystemODE(a, config); 
    // Получаем целевые значения Y*
    TargetY ty = calculateTargetStates(a, config); 
    
    double f = 0.0;
    // Вычисляем нашу свертку J (сумма ξi * Fi(a))
    for(size_t i = 0; i < 5; i++){
        f += std::abs(y[i] - ty[i]) * KSI[i]; // Ошибки для Y1, Y2, Y3, Y4, Y6
    }
    // Ошибка для Y8 (интегрального критерия)
    f += std::abs(y[5]) * KSI[5];
    return f;
}