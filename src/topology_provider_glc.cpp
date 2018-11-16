/*
 * topology_provider_glc.cpp
 *
 *  Created on: Oct 16, 2014
 *      Author: mazuran
 */

#include "topology_provider_glc.h"
#include "pseudo_chow_liu.h"
#include <g2o/types/slam2d/vertex_se2.h>
#include <g2o/types/slam3d/vertex_se3.h>
#include "glc_edge.h"
#include "glc_reparam_se2.h"
#include "glc_reparam_se3.h"
#include <Eigen/Eigen>
#include <Eigen/LU>

static const double glc_eps = 1e-8;

TopologyProviderGLC::TopologyProviderGLC() : TopologyProviderBase()
{
    _okvertices.insert(typeinfo.vertex(typeid(g2o::VertexSE2)));
    _okvertices.insert(typeinfo.vertex(typeid(g2o::VertexSE3)));
}

bool TopologyProviderGLC::applicable(
        const VertexInfoSet &vinfo,
        const EdgeInfoSet &einfo) {
    return setDifference(vinfo, _okvertices).size() == 0;
}

static GLCReparam *newReparam(bool se3)
{
    if(se3) {
        return new GLCReparamSE3;
    } else {
        return new GLCReparamSE2ISAM;
    }
}

/* TODO: shamelessly compied from isam */
Eigen::MatrixXd posdef_pinv(
        const Eigen::MatrixXd &a,
        double eps = std::numeric_limits<double>::epsilon()) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(a);
    if (eig.info() != Eigen::Success) {
        std::cout<<"mat:\n"<<a<<"\n";
        return Eigen::MatrixXd();
    }

    double tolerance = eps * std::max(a.cols(), a.rows()) * eig.eigenvalues().array().abs().maxCoeff();

    Eigen::MatrixXd mat = eig.eigenvectors() * (eig.eigenvalues().array() > tolerance).select(eig.eigenvalues().array().inverse(), 0).matrix().asDiagonal() * eig.eigenvectors().transpose();

    return mat;
}


static g2o::MatrixXD glc_chol(const g2o::MatrixXD &J, const g2o::MatrixXD &m)
{
    double eps = glc_eps * m.rows();
    int i = 0;
    Eigen::PartialPivLU<g2o::MatrixXD> lu(J);
    g2o::MatrixXD invJ = lu.solve(g2o::MatrixXD::Identity(J.rows(), J.cols()));
    g2o::MatrixXD m2 = invJ.transpose() * m * invJ;
    Eigen::SelfAdjointEigenSolver<g2o::MatrixXD> eig(m2);
    while(i < eig.eigenvalues().rows() && eig.eigenvalues()[i] < glc_eps) i++;

    return eig.eigenvectors().rightCols(m.cols() - i) *
            eig.eigenvalues().tail(m.cols() - i).array().sqrt().matrix().asDiagonal();
}

GLCEdge *TopologyProviderGLC::getEdge(
        const Eigen::MatrixXd &targetInfo,
        const g2o::OptimizableGraph::VertexContainer &vertices) const
{
    g2o::HyperGraph::VertexContainer vc(vertices.begin(), vertices.end());
    bool se3 = dynamic_cast<g2o::VertexSE3 *>(vertices.front());
    GLCEdge *edge = new GLCEdge;
    GLCReparam *reparam = newReparam(se3);
    g2o::VectorXD meas = reparam->reparametrize(vc);
    g2o::MatrixXD J = reparam->jacobian(vc, meas);
    g2o::MatrixXD W = glc_chol(J, targetInfo).transpose();

    if(W.rows() == 0) {
        delete reparam;
        delete edge;
        return NULL;
    }

    edge->setReparam(reparam);
    edge->setDimension(W.rows(), W.cols());
    edge->setLinearWeight(W);
    edge->vertices() = vc;
    edge->computeMeasurement();
    edge->information() = g2o::MatrixXD::Identity(W.rows(), W.rows());
    return edge;
}

g2o::OptimizableGraph::EdgeContainer TopologyProviderGLC::topology(
        const Eigen::MatrixXd &information,
        const g2o::OptimizableGraph::VertexContainer &vertices)
{
    g2o::OptimizableGraph::EdgeContainer newedges;
    assert(vertices.size() > 0 && "empty vertex container");

    assert((_opts.topology == SparsityOptions::Dense ||
            _opts.topology == SparsityOptions::Tree) &&
            "only dense or tree topology for GLC");
    assert(_opts.linPoint == SparsityOptions::Global &&
            "only global linearization point for GLC");

    if(vertices.size() == 1) {
        /* Trivial case */
        GLCEdge *edge = getEdge(information, vertices);
        if(edge) {
//            assert(false && "Should not happen in my tests, remove for final release");
            newedges.push_back(edge);
        }
    } else if(_opts.topology == SparsityOptions::Dense) {
        GLCEdge *edge = getEdge(information, vertices);
        if(edge) {
            newedges.push_back(edge);
        }
    } else /* if(_opts.topology == SparsityOptions::Tree) */ {
        PseudoChowLiu cl;
        cl.setSparsityOptions(_opts);
        cl.setTargetInformation(information);
        cl.setVertices(vertices);
        cl.computeSparsityPattern();
        PseudoChowLiu::SparsityPattern sp = cl.getSparsityPattern();

        /* Root node */
        Eigen::MatrixXd rootInfo = cl.marginal(sp.front().front().first);
        GLCEdge *edge = getEdge(rootInfo, g2o::OptimizableGraph::VertexContainer(
                1, sp.front().front().first));
        if(edge) {
//            assert(false && "Should not happen in my tests, remove for final release");
            newedges.push_back(edge);
        }

        for(PseudoChowLiu::CorrelatedSkeletonTree &tree: sp) {
            assert(tree.size() == 1 && "not a simple tree");
            int d1 = static_cast<g2o::OptimizableGraph::Vertex *>(tree.front().first)->dimension();
            int d2 = static_cast<g2o::OptimizableGraph::Vertex *>(tree.front().second)->dimension();
            Eigen::MatrixXd jointInfo = cl.jointMarginal(tree.front().first, tree.front().second);
            Eigen::MatrixXd target(d1 + d2, d1 + d2);

//            std::cout<<"joint info:\n";
//            std::cout<<"row: "<<jointInfo.rows()<<", col: "<<jointInfo.cols()<<"\n";
//            std::cout<<jointInfo<<"\n";

            /**
             * Canyu Le
             * WARNING: jointInfo may be inf matrix due to numerical problems, which will lead to crash since posdef_pinv() will return a empty matrix.
             * */
            // -----------------------
//            std::cout<< jointInfo<<"\n";
            auto b1 = jointInfo.block(d1, 0, d2, d1);
            auto b2 = posdef_pinv(jointInfo.block(0, 0, d1, d1));
            auto b3 = jointInfo.block(0, d1, d1, d2);
//            std::cout<<"b1: "<<b1.rows()<<", "<<b1.cols()<<"\n";
//            std::cout<<"b2: "<<b2.rows()<<", "<<b2.cols()<<"\n";
//            std::cout<<"b3: "<<b3.rows()<<", "<<b3.cols()<<"\n";
//            std::cout<<"\n";
            Eigen::MatrixXd m4;
            m4 = b1 * b2 * b3;
            target <<
                    jointInfo.block(0, 0, d1, d1), jointInfo.block(0, d1, d1, d2),
                    jointInfo.block(d1, 0, d2, d1), m4;
            // -----------------------

            g2o::OptimizableGraph::VertexContainer vc;
            vc.push_back(tree.front().first);
            vc.push_back(tree.front().second);
            GLCEdge *edge = getEdge(target.selfadjointView<Eigen::Upper>(), vc);
            if(edge) {
                newedges.push_back(edge);
            } else {
//                assert(false && "Should not happen in my tests, remove for final release");
            }
        }
    }
    return newedges;
}


