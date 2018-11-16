/*
 * se2_compatibility.h
 *
 *  Created on: Jan 9, 2014
 *      Author: mazuran
 */

#ifndef SE2_COMPATIBILITY_H_
#define SE2_COMPATIBILITY_H_

#include <g2o/types/slam2d/edge_se2.h>

typedef g2o::VertexSE2 VertexSE2ISAM;

/*
 * The implementation of g2o::EdgeSE2 computes the edge error on the SE(2) manifold, while isam
 * does it in the Euclidean space (with normalized angle). EdgeSE2ISAM subclasses g2o::EdgeSE2
 * in order to take this into account and provide a fair comparison between g2o and isam.
 */
class EdgeSE2ISAM : public g2o::EdgeSE2
{
public:
    EdgeSE2ISAM() : g2o::EdgeSE2() {}
    virtual ~EdgeSE2ISAM() {}

    void computeError()
    {
        const VertexSE2ISAM *v1 = static_cast<const VertexSE2ISAM *>(_vertices[0]);
        const VertexSE2ISAM *v2 = static_cast<const VertexSE2ISAM *>(_vertices[1]);
        g2o::SE2 delta = v1->estimate().inverse() * v2->estimate();
        _error = delta.toVector() - measurement().toVector();
        _error[2] = g2o::normalize_theta(_error[2]);
    }

    void linearizeOplus()
    {
        const VertexSE2ISAM *vi = static_cast<const VertexSE2ISAM *>(_vertices[0]);
        const VertexSE2ISAM *vj = static_cast<const VertexSE2ISAM *>(_vertices[1]);
        double thetai = vi->estimate().rotation().angle();

        Eigen::Vector2d dt = vj->estimate().translation() - vi->estimate().translation();
        double si = sin(thetai), ci = cos(thetai);

        _jacobianOplusXi(0, 0) = -ci; _jacobianOplusXi(0, 1) = -si; _jacobianOplusXi(0, 2) = -si*dt.x()+ci*dt.y();
        _jacobianOplusXi(1, 0) =  si; _jacobianOplusXi(1, 1) = -ci; _jacobianOplusXi(1, 2) = -ci*dt.x()-si*dt.y();
        _jacobianOplusXi(2, 0) =  0;  _jacobianOplusXi(2, 1) = 0;   _jacobianOplusXi(2, 2) = -1;

        _jacobianOplusXj(0, 0) =  ci; _jacobianOplusXj(0, 1) = si;  _jacobianOplusXj(0, 2) = 0;
        _jacobianOplusXj(1, 0) = -si; _jacobianOplusXj(1, 1) = ci;  _jacobianOplusXj(1, 2) = 0;
        _jacobianOplusXj(2, 0) =  0;  _jacobianOplusXj(2, 1) = 0;   _jacobianOplusXj(2, 2) = 1;
    }
};

#endif /* SE2_COMPATIBILITY_H_ */
