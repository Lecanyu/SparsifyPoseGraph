/*
 * glc_edge.h
 *
 *  Created on: 14/ott/2014
 *      Author: Mladen Mazuran
 */

#ifndef GLC_EDGE_H_
#define GLC_EDGE_H_

#include <g2o/core/base_multi_edge.h>
#include "glc_reparam.h"

class GLCEdge : public g2o::BaseMultiEdge<-1, g2o::VectorXD> {
public:
    typedef g2o::BaseMultiEdge<-1, g2o::VectorXD> Base;

    GLCEdge();
    virtual ~GLCEdge();

    void computeMeasurement();

    void setReparam(GLCReparam *reparam) {
        delete _reparam;
        _reparam = reparam;
    }

    const GLCReparam *reparam() const {
        return _reparam;
    }

    void setDimension(int errsize, int meassize) {
        resize(meassize);
        _dimension = errsize;
        _information.resize(errsize, errsize);
        _error.resize(errsize, 1);
        _measurement.resize(meassize, 1);
        _W.resize(errsize, meassize);
    }

    const g2o::MatrixXD &linearWeight() const {
        return _W;
    }

    void setLinearWeight(const g2o::MatrixXD &W) {
        _W = W;
    }

    virtual void computeError();

    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;

    virtual bool setMeasurementFromState();

    virtual void initialEstimate(
            const g2o::HyperGraph::VertexSet &fixed,
            g2o::OptimizableGraph::Vertex *toEstimate);
    virtual double initialEstimatePossible(
            const g2o::HyperGraph::VertexSet &fixed,
            g2o::OptimizableGraph::Vertex *toEstimate);

    virtual void linearizeOplus();

private:
    GLCReparam *_reparam;
    g2o::MatrixXD _W;
};

#endif /* GLC_EDGE_H_ */
