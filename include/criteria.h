#ifndef CRITERIA_H
#define CRITERIA_H
#include <vector>
#include "dynamics.h"

// Константы весов (ξ) из файла Постановки (ξ1...ξ6)
const std::vector<double> KSI = {0.007015, 0.001637, 0.008185, 0.000982, 0.491090, 0.491090};

double calculateF(const std::vector<double>& a, const AeroConstants& config);
#endif