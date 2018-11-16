/*
 * topology_provider_glc.h
 *
 *  Created on: Oct 16, 2014
 *      Author: mazuran
 */

#ifndef TOPOLOGY_PROVIDER_GLC_H_
#define TOPOLOGY_PROVIDER_GLC_H_

#include "topology_provider_base.h"

class GLCEdge;

class TopologyProviderGLC : public TopologyProviderBase
{
public:

    TopologyProviderGLC();
    virtual ~TopologyProviderGLC() {}

    virtual bool applicable(
            const VertexInfoSet &vinfo,
            const EdgeInfoSet &einfo);

    virtual g2o::OptimizableGraph::EdgeContainer topology(
            const Eigen::MatrixXd &information,
            const g2o::OptimizableGraph::VertexContainer &vertices);

    virtual bool requiresOptimization() const { return false; }

    GLCEdge *getEdge(
            const Eigen::MatrixXd &targetInfo,
            const g2o::OptimizableGraph::VertexContainer &vertices) const;
};


#endif /* TOPOLOGY_PROVIDER_GLC_H_ */
