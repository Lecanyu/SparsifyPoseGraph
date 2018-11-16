/*
 * topology_provider_binary.hpp
 *
 *  Created on: Oct 16, 2014
 *      Author: mazuran
 */

#ifndef TOPOLOGY_PROVIDER_BINARY_HPP_
#define TOPOLOGY_PROVIDER_BINARY_HPP_

#include "multi_edge_correlated.h"
#include "pseudo_chow_liu.h"

template <typename V, typename E>
TopologyProviderBinary<V, E>::TopologyProviderBinary() : TopologyProviderBase()
{
    _okvertices.insert(typeinfo.vertex(typeid(V)));
    _okedges.insert(typeinfo.edge(typeid(E)));
    _okedges.insert(typeinfo.edge(typeid(MultiEdgeCorrelated<E>)));
}

template <typename V, typename E>
g2o::OptimizableGraph::EdgeContainer TopologyProviderBinary<V, E>::topology(
        const Eigen::MatrixXd &information,
        const g2o::OptimizableGraph::VertexContainer &vertices)
{
    g2o::OptimizableGraph::EdgeContainer newedges;
    if(vertices.size() < 2) return newedges;
    PseudoChowLiu cl;
    cl.setSparsityOptions(_opts);
    cl.setTargetInformation(information);
    cl.setVertices(vertices);
    cl.computeSparsityPattern();
    PseudoChowLiu::SparsityPattern sp = cl.getSparsityPattern();
    newedges.reserve(sp.size());
    for(PseudoChowLiu::SparsityPattern::const_iterator it = sp.begin(); it != sp.end(); ++it) {
        if(it->size() == 1) {
            V *v1 = dynamic_cast<V *>(it->front().first);
            V *v2 = dynamic_cast<V *>(it->front().second);

            assert(v1 && v2 && "Unexpected vertex type");

            E *e = new E;
            e->setVertex(0, v1);
            e->setVertex(1, v2);
            e->setMeasurementFromState();
            newedges.push_back(e);
        } else {
            E *dummy = new E;
            MultiEdgeCorrelated<E> *e = new MultiEdgeCorrelated<E>;
            e->setMeasurementCount(it->size());

            for(PseudoChowLiu::CorrelatedSkeletonTree::const_iterator its =
                    it->begin(); its != it->end(); ++its) {
                V *v1 = dynamic_cast<V *>(its->first);
                V *v2 = dynamic_cast<V *>(its->second);

                assert(v1 && v2 && "Unexpected vertex type");

                dummy->setVertex(0, v1);
                dummy->setVertex(1, v2);
                dummy->setMeasurementFromState();
                e->addMeasurement(g2o::HyperGraph::VertexContainer({v1, v2}), dummy->measurement());
            }
            newedges.push_back(e);
            delete dummy;
        }
    }
    return newedges;
}


#endif /* TOPOLOGY_PROVIDER_BINARY_HPP_ */

