/*
 * test_multi_edge.cpp
 *
 *  Created on: Oct 15, 2014
 *      Author: mazuran
 */

#include <g2o/types/slam3d/edge_se3.h>
#include <g2o/types/slam2d/edge_se2.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/robust_kernel_impl.h>
#include <g2o/core/factory.h>
#include "info_block_solver.h"
#include "multi_edge_correlated.h"
#include "utils.h"

typedef MultiEdgeCorrelated<g2o::EdgeSE2> MultiSE2;
typedef MultiEdgeCorrelated<g2o::EdgeSE3> MultiSE3;

typedef InfoBlockSolver<g2o::BlockSolverTraits<-1, -1> >  SlamBlockSolver;
typedef g2o::LinearSolverCholmod<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;
typedef g2o::OptimizationAlgorithmLevenberg SlamOptimizationAlgorithm;

G2O_REGISTER_TYPE(EDGE_MULTI_SE2, MultiSE2);

int main(int argc, char *argv[])
{
    g2o::SparseOptimizer *so = new g2o::SparseOptimizer;
    SlamLinearSolver *linearSolver = new SlamLinearSolver;
    linearSolver->setBlockOrdering(false);
    SlamBlockSolver *blockSolver = new SlamBlockSolver(linearSolver);

    so->setAlgorithm(new SlamOptimizationAlgorithm(blockSolver));
    so->setVerbose(true);

    g2o::VertexSE2 *v1 = new g2o::VertexSE2;
    g2o::VertexSE2 *v2 = new g2o::VertexSE2;
    g2o::VertexSE2 *v3 = new g2o::VertexSE2;
    v1->setId(0);
    v2->setId(1);
    v3->setId(2);
    v1->setEstimate(g2o::SE2(Eigen::Vector3d::Random()));
    v2->setEstimate(g2o::SE2(Eigen::Vector3d::Random()));
    v3->setEstimate(g2o::SE2(Eigen::Vector3d::Random()));

    MultiSE2 *mse2 = new MultiSE2;
    mse2->setMeasurementCount(2);
    mse2->addMeasurement({v1, v2}, g2o::SE2(Eigen::Vector3d::Random()));
    mse2->addMeasurement({v1, v3}, g2o::SE2(Eigen::Vector3d::Random()));

    Eigen::MatrixXd A = Eigen::MatrixXd::Random(6, 6);
    A = (A * A.transpose()).eval();
    A = 0.5 * (A + A.transpose()).eval();

    mse2->setInformation(A);

    so->addVertex(v1);
    so->addVertex(v2);
    so->addVertex(v3);
    so->addEdge(mse2);

    Eigen::VectorXd v(9);
    v.segment<3>(0) = v1->estimate().toVector();
    v.segment<3>(3) = v2->estimate().toVector();
    v.segment<3>(6) = v3->estimate().toVector();
    std::cout << "before: " << v.transpose() << std::endl;

    v1->setFixed(true);
    so->initializeOptimization();
    so->optimize(10);

    so->save("test_multi_edge.g2o");

    v.segment<3>(0) = v1->estimate().toVector();
    v.segment<3>(3) = v2->estimate().toVector();
    v.segment<3>(6) = v3->estimate().toVector();

    std::cout << "after: " << v.transpose() << std::endl << std::endl;
    std::cout << "m1 = " << mse2->measurement().front().toVector().transpose() << std::endl;
    std::cout << "m2 = " << mse2->measurement().back().toVector().transpose() << std::endl;
    std::cout << "minfo = " << std::endl << mse2->information() << std::endl;

    so->removeEdge(mse2);
    so->removeVertex(v1);
    so->removeVertex(v2);
    so->removeVertex(v3);
    so->load("test_multi_edge.g2o");

    std::cout << std::endl;
    std::cout << "m1 = " << static_cast<MultiSE2 *>(*so->edges().begin())->measurement().front().toVector().transpose() << std::endl;
    std::cout << "m2 = " << static_cast<MultiSE2 *>(*so->edges().begin())->measurement().back().toVector().transpose() << std::endl;
    std::cout << "minfo = " << std::endl << static_cast<MultiSE2 *>(*so->edges().begin())->information() << std::endl;

    delete so;

    return 0;
}
