/*
 * pqn_optimizer_options.h
 *
 *  Created on: Feb 14, 2014
 *      Author: mazuran
 */

#ifndef PQN_OPTIMIZER_OPTIONS_H_
#define PQN_OPTIMIZER_OPTIONS_H_

struct PQNOptimizerOptions {
    bool verbose;
    bool usePreconditioner;
    bool useHessian;
    int maxIters;
    int lbfgsCorrections;
    double sufficientDecrease;
    double curvatureCondition;
    double optimalityTolerance;

    PQNOptimizerOptions() :
        verbose(false),
        usePreconditioner(false),
        useHessian(false),
        maxIters(0),
        lbfgsCorrections(10),
        sufficientDecrease(1e-6),
        curvatureCondition(0.9),
        optimalityTolerance(1e-12) {}
};

#endif /* PQN_OPTIMIZER_OPTIONS_H_ */
