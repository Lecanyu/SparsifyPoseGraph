/*
 * topology_provider_se3.h
 *
 *  Created on: Oct 16, 2014
 *      Author: mazuran
 */

#ifndef TOPOLOGY_PROVIDER_SE3_H_
#define TOPOLOGY_PROVIDER_SE3_H_

#include "topology_provider_binary.h"
#include <g2o/types/slam3d/vertex_se3.h>
#include <g2o/types/slam3d/edge_se3.h>
#include "se3_compatibility.h"

typedef TopologyProviderBinary<g2o::VertexSE3, g2o::EdgeSE3> TopologyProviderSE3;
typedef TopologyProviderBinary<g2o::VertexSE3, EdgeSE3ISAM> TopologyProviderSE3ISAM;

#endif /* TOPOLOGY_PROVIDER_SE3_H_ */
