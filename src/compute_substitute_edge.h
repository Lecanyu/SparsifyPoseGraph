/*
 * compute_substitute_edge.h
 *
 *  Created on: Oct 23, 2014
 *      Author: mazuran
 */

#include "graph_wrapper.h"
#include <set>
#include <Eigen/Core>

void computeSubstituteEdge(
        GraphWrapper *gw, const std::set<int> &marginalized, int maxid,
        int &from, int &to, IsometryXd &edgemeas, Eigen::MatrixXd &edgeinfo);


