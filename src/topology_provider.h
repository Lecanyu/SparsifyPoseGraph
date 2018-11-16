/*
 * topology_provider.h
 *
 *  Created on: Oct 16, 2014
 *      Author: mazuran
 */

#ifndef TOPOLOGY_PROVIDER_H_
#define TOPOLOGY_PROVIDER_H_

#include <set>
#include <Eigen/Core>
#include "type_info.h"
#include "sparsity_options.h"

class TopologyProvider
{
public:
    typedef std::set<const VertexInfo *> VertexInfoSet;
    typedef std::set<const EdgeInfo *> EdgeInfoSet;

    virtual ~TopologyProvider() {}

    virtual bool applicable(
            const VertexInfoSet &vinfo,
            const EdgeInfoSet &einfo) = 0;

    virtual g2o::OptimizableGraph::EdgeContainer topology(
            const Eigen::MatrixXd &information,
            const g2o::OptimizableGraph::VertexContainer &vertices) = 0;

    virtual void setSparsityOptions(const SparsityOptions &opts) { _opts = opts; }
    virtual bool requiresOptimization() const = 0;

protected:
    SparsityOptions _opts;
};



#endif /* TOPOLOGY_PROVIDER_H_ */
