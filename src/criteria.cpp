#include "criteria.h"
#include "target_states.h"
#include "stability.h"
#include "dynamics.h"
#include <cmath>

double calculateF(const std::vector<double>& a, const AeroConstants& config) {
    SimulationResult y = solveSystemODE(a, config); // 
    TargetY ty = calculateTargetStates(a, config); //
    double f = 0.0;
    for(size_t i = 0; i < 5; i++){
        f += abs(y[i] - ty[i]) * KSI[i]; // F1..F5
    }
    // F6
    
    f += abs(y[5]) * KSI[5];
    return f;
}
