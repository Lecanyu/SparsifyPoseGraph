/*
 * compute_substitute_edge.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: mazuran
 */

#include "compute_substitute_edge.h"
#include <deque>
#include <limits>
#include "isometryxd.h"

void computeSubstituteEdge(
        GraphWrapper *gw, const std::set<int> &marginalized, int maxid,
        int &from, int &to, IsometryXd &edgemeas, Eigen::MatrixXd &edgeinfo)
{
    std::set<int> visited;
    std::deque<std::set<int> > frontiers;
    std::set<int> newFrontier;
    int minid = std::numeric_limits<int>::max();

    int toConnect = std::max(from, to);
    int toReplace = std::min(from, to);

    newFrontier.insert(toReplace);
    visited.insert(toConnect);
    visited.insert(toReplace);
    do {
        frontiers.push_back(newFrontier);
        newFrontier.clear();
        for(int r: frontiers.back()) {
            if(marginalized.count(r) == 0 && r != from && r != to) {
                minid = std::min(minid, r);
            } else {
                auto edges = gw->vertex(r)->edges();
                visited.insert(r);
                for(auto e: edges) {
                    auto vertices = e->vertices();
                    if(vertices.size() == 2) {
                        int idother = vertices[0]->id() == r ? vertices[1]->id() : vertices[0]->id();
                        if(visited.count(idother) == 0 && idother <= maxid && idother != 0) {
                            newFrontier.insert(idother);
                        }
                    }
                }
            }
        }
    } while(minid == std::numeric_limits<int>::max());
    visited.clear(); visited.insert(toConnect); /* Reuse visited var, don't need it anymore */
    frontiers.push_front(visited);
    frontiers.pop_back();

    int dim = gw->vertex(toConnect)->is2d() ? 3 : 6;
    Eigen::MatrixXd covsum = Eigen::MatrixXd::Zero(dim, dim);
    IsometryXd meas(gw->vertex(toConnect)->is2d());
    int reach = minid;

    while(!frontiers.empty()) {
        std::set<int> lastFrontier = frontiers.back();
        frontiers.pop_back();
        auto edges = gw->vertex(reach)->edges();
        for(auto e: edges) {
            auto endpoints = e->vertices();
            if(endpoints.size() != 2) continue;
            if(lastFrontier.count(endpoints[0]->id()) || lastFrontier.count(endpoints[1]->id())) {
                covsum += e->information().inverse();
                if(from == toConnect) {
                    if(endpoints[1]->id() == reach) {
                        meas = e->measurement() * meas;
                    } else {
                        meas = e->measurement().inverse() * meas;
                    }
                } else {
                    if(endpoints[1]->id() == reach) {
                        meas *= e->measurement().inverse();
                    } else {
                        meas *= e->measurement();
                    }
                }
                reach = endpoints[1]->id() == reach ? endpoints[0]->id() : endpoints[1]->id();
                break;
            }
        }
    }


    edgeinfo = covsum.inverse();
    edgeinfo = 0.5 * (edgeinfo + edgeinfo.transpose()).eval();
    edgemeas = meas;

    if(from == toConnect) {
        to = minid;
    } else {
        from = minid;
    }
}


