/*
 * pqn_general.cpp
 *
 *  Created on: Feb 14, 2014
 *      Author: mazuran
 */

#include "pqn_general.h"
#include <cmath>

bool invalid(double v)
{
    return std::isnan(v) || std::isinf(v);
}

bool invalid(const Eigen::VectorXd &v)
{
    for(int i = 0; i < v.rows(); i++) {
        if(invalid(v[i])) return true;
    }
    return false;
}

double twoPointCubicInterpolation(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2)
{
    double d1 = p1[2] + p2[2] - 3 * (p1[1] - p2[1]) / (p1[0] - p2[0]);
    double diff = d1 * d1 - p1[2] * p2[2];
    if(diff >= 0) {
        double d2 = std::sqrt(diff);
        double t = p2[0] - (p2[0] - p1[0]) * (p2[2] + d2 - d1) / (p2[2] - p1[2] + 2 * d2);
        return std::min(std::max(t, p1[0]), p2[0]);
    } else {
        return 0.5 * (p1[0] + p2[0]);
    }
}
