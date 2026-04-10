#ifndef TARGET_STATES_H
#define TARGET_STATES_H
#include <vector>
#include "dynamics.h"

using TargetY = std::vector<double>;
TargetY calculateTargetStates(const std::vector<double>& a, const AeroConstants& config);
#endif