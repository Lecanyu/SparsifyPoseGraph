/*
 * test_logdet.cpp
 *
 *  Created on: 18/ott/2014
 *      Author: Mladen Mazuran
 */

#include <Eigen/Core>
#include <iostream>
#include "logdet_function.h"

int main(int argc, char *argv[]) {
    Eigen::MatrixXd J11, J12, J21, J22, T;

    J11 = Eigen::MatrixXd::Random(3, 3);
    J12 = Eigen::MatrixXd::Random(3, 3);
    J21 = Eigen::MatrixXd::Random(3, 3);
    J22 = Eigen::MatrixXd::Random(3, 3);
    T   = Eigen::MatrixXd::Random(6, 3);
    T   = (T * T.transpose()).eval();
    T   = 0.5 * (T + T.transpose()).eval();

    JacobianEntry e11, e12, e21, e22;
    e11.first = J11;
    e12.first = J12;
    e21.first = J21;
    e22.first = J22;

    e11.second = 0;
    e12.second = 3;
    e21.second = 0;
    e22.second = 3;

    MeasurementJacobian mj1, mj2;
    mj1.push_back(e11);
    mj1.push_back(e12);
    mj2.push_back(e21);
    mj2.push_back(e22);

    JacobianMapping jm;
    jm.push_back(mj1);
    jm.push_back(mj2);

    LogdetFunctionWithConstraints fun(jm, T);
//    Eigen::VectorXd guess(12), g;
//    guess <<
//            0.3, 0.1, 0, 0.3, 0.1, 0.3,
//            0.3, 0.1, 0, 0.3, 0.1, 0.3;
    Eigen::VectorXd guess(18), g;
    guess <<
            0.3, 0.1, 0.0, 0.1, 0.3, 0.1,
            0.0, 0.1, 0.3, 0.3, 0.1, 0.0,
            0.1, 0.3, 0.1, 0.0, 0.1, 0.3;

    fun.setRho(1);

    std::cout << J11 << std::endl << std::endl;
    std::cout << J12 << std::endl << std::endl;
    std::cout << J21 << std::endl << std::endl;
    std::cout << J22 << std::endl << std::endl;
    std::cout << T << std::endl << std::endl;
    std::cout << fun.value(guess) << std::endl;
    fun.gradient(guess, g);
    std::list<Eigen::MatrixXd> l = fun.decondense(g);
    for(auto &m: l) {
        std::cout << std::endl << m << std::endl;
    }

    Eigen::MatrixXd H(18, 18);
    fun.hessian(guess, H);

    std::cout << std::endl << H << std::endl;

    return 0;
}


