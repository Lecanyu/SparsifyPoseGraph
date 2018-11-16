/*
 * topology_provider_se2.h
 *
 *  Created on: Oct 16, 2014
 *      Author: mazuran
 */

#ifndef TOPOLOGY_PROVIDER_SE2_H_
#define TOPOLOGY_PROVIDER_SE2_H_

#include "topology_provider_binary.h"
#include <g2o/types/slam2d/vertex_se2.h>
#include <g2o/types/slam2d/edge_se2.h>
#include "se2_compatibility.h"

typedef TopologyProviderBinary<g2o::VertexSE2, g2o::EdgeSE2> TopologyProviderSE2;
typedef TopologyProviderBinary<g2o::VertexSE2, EdgeSE2ISAM> TopologyProviderSE2ISAM;

#endif /* TOPOLOGY_PROVIDER_SE2_H_ */
