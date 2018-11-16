/*
 * edge_types.cpp
 *
 *  Created on: 14/oct/2014
 *      Author: Mladen Mazuran
 */

#include "type_info.h"
#include <g2o/types/slam2d/edge_se2.h>
#include <g2o/types/slam3d/edge_se3.h>
#include <g2o/core/factory.h>
#include "multi_edge_correlated.h"
#include "se2_compatibility.h"
#include "se3_compatibility.h"
#include "glc_edge.h"
#include "glc_reparam_se2.h"
#include "glc_reparam_se3.h"

template <typename EdgeType>
class MultiEdgeInfo : public EdgeInfoBase<MultiEdgeCorrelated<EdgeType> >
{
public:
    virtual g2o::HyperGraph::Edge *clone(
            const MultiEdgeCorrelated<EdgeType> *e) const {
        MultiEdgeCorrelated<EdgeType> *c = new MultiEdgeCorrelated<EdgeType>;
        c->setMeasurementCount(e->measurementCount());
        for(typename MultiEdgeCorrelated<EdgeType>::const_iterator it = e->begin();
                it != e->end(); ++it) {
            c->addMeasurement(it.vertices(), it.measurement());
        }
        c->setInformation(e->information());
        return c;
    }
};

class GLCEdgeInfo : public EdgeInfoBase<GLCEdge>
{
public:
    virtual g2o::HyperGraph::Edge *clone(
            const GLCEdge *e) const {
        GLCEdge *c = new GLCEdge;
        std::string tag = g2o::Factory::instance()->tag(e->reparam());
        c->setReparam(static_cast<GLCReparam *>(g2o::Factory::instance()->construct(tag)));
        c->setDimension(e->linearWeight().rows(), e->linearWeight().cols());
        c->vertices() = e->vertices();
        c->setMeasurement(e->measurement());
        c->setInformation(e->information());
        c->setLinearWeight(e->linearWeight());
        return c;
    }
};

typedef EdgeInfoBase<EdgeSE2ISAM> EdgeSE2ISAMInfo;
typedef EdgeInfoBase<EdgeSE3ISAM> EdgeSE3ISAMInfo;
typedef EdgeInfoBase<g2o::EdgeSE2> EdgeSE2Info;
typedef EdgeInfoBase<g2o::EdgeSE3> EdgeSE3Info;

typedef MultiEdgeCorrelated<EdgeSE2ISAM> MultiEdgeSE2ISAM;
typedef MultiEdgeCorrelated<EdgeSE3ISAM> MultiEdgeSE3ISAM;
typedef MultiEdgeCorrelated<g2o::EdgeSE2> MultiEdgeSE2;
typedef MultiEdgeCorrelated<g2o::EdgeSE3> MultiEdgeSE3;

typedef MultiEdgeInfo<EdgeSE2ISAM> MultiEdgeSE2ISAMInfo;
typedef MultiEdgeInfo<EdgeSE3ISAM> MultiEdgeSE3ISAMInfo;
typedef MultiEdgeInfo<g2o::EdgeSE2> MultiEdgeSE2Info;
typedef MultiEdgeInfo<g2o::EdgeSE3> MultiEdgeSE3Info;

G2O_REGISTER_TYPE(GLC_REPARAM_SE2_ISAM, GLCReparamSE2ISAM);
G2O_REGISTER_TYPE(GLC_REPARAM_SE2, GLCReparamSE2);
G2O_REGISTER_TYPE(GLC_REPARAM_SE3, GLCReparamSE3);
G2O_REGISTER_TYPE(GLC_EDGE, GLCEdge);
G2O_REGISTER_TYPE(EDGE_SE2_ISAM, EdgeSE2ISAM);
G2O_REGISTER_TYPE(MULTI_EDGE_SE2_ISAM, MultiEdgeSE2ISAM);
G2O_REGISTER_TYPE(MULTI_EDGE_SE3_ISAM, MultiEdgeSE3ISAM);
G2O_REGISTER_TYPE(MULTI_EDGE_SE2, MultiEdgeSE2);
G2O_REGISTER_TYPE(MULTI_EDGE_SE3, MultiEdgeSE3);

G2S_REGISTER_TYPE(GLCEdgeInfo);
G2S_REGISTER_TYPE(EdgeSE2ISAMInfo);
G2S_REGISTER_TYPE(EdgeSE3ISAMInfo);
G2S_REGISTER_TYPE(EdgeSE2Info);
G2S_REGISTER_TYPE(EdgeSE3Info);

G2S_REGISTER_TYPE(MultiEdgeSE2ISAMInfo);
G2S_REGISTER_TYPE(MultiEdgeSE3ISAMInfo);
G2S_REGISTER_TYPE(MultiEdgeSE2Info);
G2S_REGISTER_TYPE(MultiEdgeSE3Info);
