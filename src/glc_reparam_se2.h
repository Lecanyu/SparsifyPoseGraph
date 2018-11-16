/*
 * glc_reparam_se2.h
 *
 *  Created on: Oct 20, 2014
 *      Author: mazuran
 */

#ifndef GLC_REPARAM_SE2_H_
#define GLC_REPARAM_SE2_H_

#include "glc_reparam_binary.h"
#include <g2o/types/slam2d/edge_se2.h>
#include <g2o/types/slam2d/vertex_se2.h>
#include "se2_compatibility.h"

class GLCReparamSE2 : public GLCReparamBinary<g2o::VertexSE2, g2o::EdgeSE2>
{
public:
    virtual ~GLCReparamSE2() {}

    virtual g2o::SE2 zero() const {
        return g2o::SE2(0, 0, 0);
    }

    virtual g2o::SE2 errorToMeasurement(const g2o::VectorXD &err) const {
        return g2o::SE2(err);
    }
};

class GLCReparamSE2ISAM : public GLCReparamBinary<g2o::VertexSE2, EdgeSE2ISAM>
{
public:
    virtual ~GLCReparamSE2ISAM() {}

    virtual g2o::SE2 zero() const {
        return g2o::SE2(0, 0, 0);
    }

    virtual g2o::SE2 errorToMeasurement(const g2o::VectorXD &err) const {
        return g2o::SE2(err);
    }
};

#endif /* GLC_REPARAM_SE2_H_ */
