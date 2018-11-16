/*
 * glc_reparam_se3.h
 *
 *  Created on: Oct 20, 2014
 *      Author: mazuran
 */

#ifndef GLC_REPARAM_SE3_H_
#define GLC_REPARAM_SE3_H_

#include "glc_reparam_binary.h"
#include <g2o/types/slam3d/edge_se3.h>
#include <g2o/types/slam3d/vertex_se3.h>
#include <g2o/types/slam3d/isometry3d_mappings.h>

class GLCReparamSE3 : public GLCReparamBinary<g2o::VertexSE3, g2o::EdgeSE3>
{
public:
    virtual ~GLCReparamSE3() {}

    virtual g2o::Isometry3D zero() const {
        return g2o::Isometry3D::Identity();
    }

    virtual g2o::Isometry3D errorToMeasurement(const g2o::VectorXD &err) const {
        return g2o::internal::fromVectorMQT(err);
    }
};

#endif /* GLC_REPARAM_SE3_H_ */
