#ifndef TARGET_STATES_H
#define TARGET_STATES_H
#include <vector>
#include "dynamics.h"

struct TargetY {
    double y1 = 0.0;
    double y2 = 0.0;
    double y3 = 0.0;
    double y4 = 0.0;
    double y6 = 0.0;
};
TargetY calculateTargetStates(const std::vector<double>& a, const AeroConstants& config);
#endif