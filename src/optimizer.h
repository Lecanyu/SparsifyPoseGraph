/*
 * optimizer.h
 *
 *  Created on: 08/dic/2013
 *      Author: Mladen Mazuran
 */

#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_

#include <vector>
#include <list>
#include <Eigen/Core>
#include <g2o/core/eigen_types.h>

/* Jacobian wrt vertex, start index of vertex in info matrix */
typedef std::pair<g2o::MatrixXD, int> JacobianEntry;
typedef std::list<JacobianEntry> MeasurementJacobian;
typedef std::list<MeasurementJacobian> JacobianMapping;

std::list<g2o::MatrixXD> optimizeInformation(
        const JacobianMapping &mapping,
        const g2o::MatrixXD &targetInformation);

#endif /* OPTIMIZER_H_ */
