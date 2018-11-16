/*
 * glc_reparam.h
 *
 *  Created on: Oct 20, 2014
 *      Author: mazuran
 */

#ifndef GLC_REPARAM_H_
#define GLC_REPARAM_H_

#include <g2o/core/hyper_graph.h>
#include <g2o/core/eigen_types.h>
#include <g2o/core/creators.h>

class GLCReparam : public g2o::HyperGraph::HyperGraphElement {
public:
    virtual ~GLCReparam() {}

    virtual g2o::VectorXD reparametrize(
            const g2o::HyperGraph::VertexContainer &vertices) = 0;

    virtual g2o::VectorXD reparametrize(
            const g2o::HyperGraph::VertexContainer &vertices,
            const g2o::VectorXD &measure) = 0;

    virtual g2o::MatrixXD jacobian(
            const g2o::HyperGraph::VertexContainer &vertices) = 0;

    virtual g2o::MatrixXD jacobian(
            const g2o::HyperGraph::VertexContainer &vertices,
            const g2o::VectorXD &meas) = 0;

    virtual g2o::HyperGraph::HyperGraphElementType elementType() const {
        return g2o::HyperGraph::HGET_DATA;
    }
};



#endif /* GLC_REPARAM_H_ */
