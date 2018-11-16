/*
 * pqn_general.h
 *
 *  Created on: Feb 14, 2014
 *      Author: mazuran
 */

#ifndef PQN_GENERAL_H_
#define PQN_GENERAL_H_

#include <Eigen/Core>

bool invalid(double v);
bool invalid(const Eigen::VectorXd &v);
double twoPointCubicInterpolation(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2);

#endif /* PQN_GENERAL_H_ */
