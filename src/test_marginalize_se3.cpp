/*
 * test_marginalize_se3.cpp
 *
 *  Created on: Oct 9, 2014
 *      Author: mazuran
 */

#include <iostream>
#include "graph_wrapper_g2o.h"

g2o::Vector6d euler(double x, double y, double z, double yaw, double pitch, double roll)
{
    g2o::Vector6d euler;
    euler << x, y, z, yaw, pitch, roll;
    return euler;
}

int main(int argc, char *argv[])
{
    GraphWrapperG2O w(false);
    g2o::Vector6d m = euler(0.1, 0.1, 0.1, 0.1, 0.1, 0.1);
    g2o::Vector6d m1 = m + 0.05 * g2o::Vector6d::Random();
    g2o::Vector6d m2 = m + 0.05 * g2o::Vector6d::Random();

    w.addVertex(0, IsometryXd(false));
    w.addVertex(1, IsometryXd(m1, IsometryXd::EulerAnglesG2O));
    w.addVertex(2, IsometryXd(m1, IsometryXd::EulerAnglesG2O) * IsometryXd(m2, IsometryXd::EulerAnglesG2O));
    w.addEdge(0, 1, IsometryXd(m, IsometryXd::EulerAnglesG2O), Eigen::MatrixXd::Identity(6, 6));
    w.addEdge(1, 2, IsometryXd(m, IsometryXd::EulerAnglesG2O), Eigen::MatrixXd::Identity(6, 6));

    w.optimize();
    std::cout << w.estimate() << std::endl << std::endl;

    std::vector<int> marg; marg.push_back(1);
    GraphWrapperG2O *mw = w.clonePortion(2);
    SparsityOptions opts;
    opts.linPoint = SparsityOptions::Local;
    opts.topology = SparsityOptions::Tree;
    mw->marginalize(marg, opts);
    Eigen::MatrixXd cov = mw->covariance();

//    mw->debugPrint(std::cout);
//    std::cout<<std::endl;
//    mw->printStats(std::cout);
//    std::cout<<std::endl;
    std::cout << mw->estimate() << std::endl << std::endl;

    std::cout << "kld = " << w.kullbackLeibler(mw) << std::endl;

    delete mw;
}


