/*
 * glc_reparam_binary.hpp
 *
 *  Created on: Oct 20, 2014
 *      Author: mazuran
 */

#ifndef GLC_REPARAM_BINARY_HPP_
#define GLC_REPARAM_BINARY_HPP_

#include <g2o/core/jacobian_workspace.h>

template <typename VertexType, typename EdgeType>
GLCReparamBinary<VertexType,EdgeType>::GLCReparamBinary() :
    GLCReparam(), _mock(new EdgeType), _zero(new VertexType)
{
    _mock->setInformation(g2o::MatrixXD(int(EdgeType::Dimension), int(EdgeType::Dimension)));
}

template <typename VertexType, typename EdgeType>
GLCReparamBinary<VertexType,EdgeType>::~GLCReparamBinary()
{
    delete _mock;
    delete _zero;
}

template <typename VertexType, typename EdgeType>
g2o::VectorXD GLCReparamBinary<VertexType,EdgeType>::reparametrize(
        const g2o::HyperGraph::VertexContainer &vertices)
{
    return reparametrize(vertices, g2o::VectorXD::Zero(vertices.size() * EdgeType::Dimension));
}

template <typename VertexType, typename EdgeType>
g2o::VectorXD GLCReparamBinary<VertexType,EdgeType>::reparametrize(
        const g2o::HyperGraph::VertexContainer &vertices,
        const g2o::VectorXD &meas)
{
    int k = EdgeType::Dimension;
    g2o::VectorXD ret(EdgeType::Dimension * vertices.size());
    assert(vertices.size() > 0 && "empty vertex container");
    VertexType *v0 = dynamic_cast<VertexType *>(vertices.front());
    assert(v0 && "wrong vertex type");

    /*
     * The GLC paper says that the first entry in the reparametrization should
     * be the inverse of the pose, but in the actual isam implementation it's
     * the pose itself. Here we replicate the isam code choice.
     */
    _zero->setEstimate(zero());
    _mock->setVertex(0, _zero);
    _mock->setVertex(1, v0);
    _mock->setMeasurement(errorToMeasurement(meas.head<EdgeType::Dimension>()));

    _mock->computeError();
    ret.head<EdgeType::Dimension>() = _mock->error();
    _mock->setVertex(0, v0);
    for(size_t i = 1; i < vertices.size(); i++, k += EdgeType::Dimension) {
        VertexType *v1 = dynamic_cast<VertexType *>(vertices[i]);
        assert(v1 && "wrong vertex type");
        _mock->setVertex(1, v1);
        _mock->setMeasurement(errorToMeasurement(meas.segment<EdgeType::Dimension>(k)));
        _mock->computeError();
        ret.segment<EdgeType::Dimension>(k) = _mock->error();
    }

    return ret;
}

template <typename VertexType, typename EdgeType>
g2o::MatrixXD GLCReparamBinary<VertexType,EdgeType>::jacobian(
            const g2o::HyperGraph::VertexContainer &vertices)
{
    return jacobian(vertices, g2o::VectorXD::Zero(vertices.size() * EdgeType::Dimension));
}

template <typename VertexType, typename EdgeType>
g2o::MatrixXD GLCReparamBinary<VertexType,EdgeType>::jacobian(
            const g2o::HyperGraph::VertexContainer &vertices,
            const g2o::VectorXD &meas)
{
    int k = EdgeType::Dimension;
    g2o::MatrixXD J = g2o::MatrixXD::Zero(
            EdgeType::Dimension * vertices.size(),
            EdgeType::Dimension * vertices.size());
    assert(vertices.size() > 0 && "empty vertex container");

    VertexType *v0 = dynamic_cast<VertexType *>(vertices.front());
    assert(v0 && "wrong vertex type");
    g2o::JacobianWorkspace jw;

    /*
     * The GLC paper says that the first entry in the reparametrization should
     * be the inverse of the pose, but in the actual isam implementation it's
     * the pose itself. Here we replicate the isam code choice.
     */
    _mock->setVertex(0, _zero);
    _mock->setVertex(1, v0);
    _mock->setMeasurement(errorToMeasurement(meas.head<EdgeType::Dimension>()));
    jw.updateSize(_mock);
    jw.allocate();

    static_cast<g2o::OptimizableGraph::Edge *>(_mock)->linearizeOplus(jw);
    J.block<EdgeType::Dimension, VertexType::Dimension>(0, 0) = _mock->jacobianOplusXj();

    _mock->setVertex(0, v0);
    for(size_t i = 1; i < vertices.size(); i++, k += EdgeType::Dimension) {
        VertexType *v1 = dynamic_cast<VertexType *>(vertices[i]);
        assert(v1 && "wrong vertex type");
        _mock->setVertex(1, v1);
        _mock->setMeasurement(errorToMeasurement(meas.segment<EdgeType::Dimension>(k)));
        _mock->linearizeOplus();
        J.block<EdgeType::Dimension, VertexType::Dimension>(k, 0) =
                        _mock->jacobianOplusXi();
        J.block<EdgeType::Dimension, VertexType::Dimension>(k, k) =
                        _mock->jacobianOplusXj();
    }

    return J;
}




#endif /* GLC_REPARAM_BINARY_HPP_ */
