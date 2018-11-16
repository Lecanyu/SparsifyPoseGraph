/*
 * optimizer.cpp
 *
 *  Created on: 08/dic/2013
 *      Author: Mladen Mazuran
 */

#include <iostream>
#include "utils.h"
#include "pqn/pqn_optimizer.h"
#include "logdet_function.h"
#include "optimizer.h"
#include <cfloat>
#include <cmath>

std::list<g2o::MatrixXD> optimizeInformation(
        const JacobianMapping &mapping,
        const g2o::MatrixXD &targetInformation)
{
    LogdetFunctionWithConstraints fun(mapping, targetInformation);
    if(fun.hasClosedFormSolution()) {
        return fun.closedFormSolution();
//        std::list<g2o::MatrixXD> sol = fun.closedFormSolution();
//        Eigen::MatrixXd infoprod = fun.informationProduct(sol);
//        double kld = fun.value(fun.condense(sol));
//        Eigen::MatrixXd diff = targetInformation - infoprod;
//        double froDiff = std::sqrt(diff.array().square().sum());
//        double maxDiff = diff.cwiseAbs().maxCoeff();
//        double maxCoeff = targetInformation.cwiseAbs().maxCoeff();
//        std::cout << std::endl << "Projected KLD = " << kld <<
//                "  Frobenius differece = " << froDiff <<
//                "  Max difference = " << maxDiff << std::endl;
//        std::cout << "Scale KLD = " << kld / maxCoeff <<
//                "  Scale Frobenius = " << froDiff / maxCoeff <<
//                "  Scale Max = " << maxDiff / maxCoeff << std::endl;
//        return sol;
    } else {
        Eigen::VectorXd x = fun.educatedGuess();

        /* Dirty interior point */
        const double startRho = 1, endRho = 5e-8, stepRho = std::sqrt(10);
        double rho;

        PQNOptimizer optimizer;
        PQNOptimizerOptions options = optimizer.options();
        /**
         * Canyu Le
         * Close the verbose flag and comment out std::cerr
         * */
        // -----------------------
        options.verbose = false;
        // -----------------------

        options.useHessian = true;
        options.optimalityTolerance = 1e-4;
        optimizer.setLineSearch(new LineSearchSimpleBacktracking);
        optimizer.setOptions(options);

        for(rho = startRho; rho >= endRho; rho /= stepRho) {
//            std::cerr << "+-----------------------+" << std::endl;
//            std::cerr << "| rho = " << std::setw(15) << rho << " |" << std::endl;
//            std::cerr << "+-----------------------+" << std::endl;
            fun.setRho(rho);
            if(rho / stepRho < endRho) {
                options.optimalityTolerance = 1e-12;
                optimizer.setOptions(options);
            }
            optimizer.optimize(fun, x);
//            std::cerr << fun.LogdetFunction::value(x) << std::endl;
        }

        double final = fun.LogdetFunction::value(x);
//        std::cerr << "Final function value: " << final << std::endl;

        if(std::isinf(final)) {
            exit(0);
        }

        return fun.decondense(x);
    }
}

