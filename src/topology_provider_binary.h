/*
 * topology_provider_binary.h
 *
 *  Created on: Oct 16, 2014
 *      Author: mazuran
 */

#ifndef TOPOLOGY_PROVIDER_BINARY_H_
#define TOPOLOGY_PROVIDER_BINARY_H_

#include "topology_provider_base.h"

template <typename V, typename E>
class TopologyProviderBinary : public TopologyProviderBase
{
public:

    TopologyProviderBinary();
    virtual ~TopologyProviderBinary() {}

    virtual g2o::OptimizableGraph::EdgeContainer topology(
            const Eigen::MatrixXd &information,
            const g2o::OptimizableGraph::VertexContainer &vertices);
};

#include "topology_provider_binary.hpp"


#endif /* TOPOLOGY_PROVIDER_BINARY_H_ */
