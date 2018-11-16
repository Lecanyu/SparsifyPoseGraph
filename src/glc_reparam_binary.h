/*
 * glc_reparam_binary.h
 *
 *  Created on: Oct 20, 2014
 *      Author: mazuran
 */

#ifndef GLC_REPARAM_BINARY_H_
#define GLC_REPARAM_BINARY_H_

#include "glc_reparam.h"
#include <list>

template <typename VertexType, typename EdgeType>
class GLCReparamBinary : public GLCReparam {
public:
    typedef typename EdgeType::Measurement Measurement;

    GLCReparamBinary();
    virtual ~GLCReparamBinary();

    virtual g2o::VectorXD reparametrize(
                const g2o::HyperGraph::VertexContainer &vertices);

    virtual g2o::VectorXD reparametrize(
                const g2o::HyperGraph::VertexContainer &vertices,
                const g2o::VectorXD &meas);

    virtual g2o::MatrixXD jacobian(
            const g2o::HyperGraph::VertexContainer &vertices);

    virtual g2o::MatrixXD jacobian(
            const g2o::HyperGraph::VertexContainer &vertices,
            const g2o::VectorXD &meas);

    virtual Measurement errorToMeasurement(
            const g2o::VectorXD &error) const = 0;

    virtual typename EdgeType::Measurement zero() const = 0;

protected:
    EdgeType *_mock;
    VertexType *_zero;
};

#include "glc_reparam_binary.hpp"

#endif /* GLC_REPARAM_BINARY_H_ */
