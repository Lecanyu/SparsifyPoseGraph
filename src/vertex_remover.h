/*
 * vertex_remover.h
 *
 *  Created on: Oct 16, 2014
 *      Author: mazuran
 */

#ifndef VERTEX_REMOVER_H_
#define VERTEX_REMOVER_H_

#include <g2o/core/sparse_optimizer.h>
#include "sparsity_options.h"
#include "topology_provider.h"
#include "optimizer.h"

#include "auto_delete_map.h"
#include "graph_wrapper.h"

class VertexRemover {
public:
    typedef std::set<g2o::OptimizableGraph::Vertex *,
            bool (*)(g2o::OptimizableGraph::Vertex *,
                    g2o::OptimizableGraph::Vertex *)> VertexSet;
    typedef g2o::HyperGraph::EdgeSet EdgeSet;
    typedef std::map<g2o::OptimizableGraph::Vertex *,
            g2o::OptimizableGraph::Vertex *> VertexMapping;
    typedef AutoDeleteMap<const g2o::HyperGraph::Edge *,   GraphWrapper::Edge *>   EdgeMap;

    VertexRemover();
    virtual ~VertexRemover();

    void setSparsityOptions(const SparsityOptions &opts);
    void setGraph(g2o::OptimizableGraph *graph) { _graph = graph; }
    void setEdgeMap(EdgeMap *edgeMap) { _edgeLookup = edgeMap; }
    void registerTopologyProvider(TopologyProvider *topology) {
        _topologies.push_back(topology);
    }

    g2o::OptimizableGraph::EdgeContainer remove(
            g2o::OptimizableGraph::Vertex *toRemove);
    g2o::OptimizableGraph::EdgeContainer remove(
            const std::vector<g2o::OptimizableGraph::Vertex *> &toRemove);

    VertexSet extendedMarkovBlanketVertices(
            g2o::OptimizableGraph::Vertex *root, const VertexSet &pickBin,
            VertexSet &picked) const;
    VertexSet markovBlanketVertices(
            g2o::OptimizableGraph::Vertex *root) const;
    EdgeSet markovBlanketEdges(
            g2o::OptimizableGraph::Vertex *root) const;

private:
    EdgeSet markovBlanketEdges(
            const VertexSet &mbVertices,
            const VertexSet &hubs) const;

    void addSubgraphVertex(
            int id, g2o::OptimizableGraph::Vertex *original,
            const g2o::OptimizableGraph::Vertex *reparam = NULL);

    void addSubgraphEdge(const g2o::HyperGraph::Edge *e);

    void buildSubgraph(
            const VertexSet &toRemove,
            const VertexSet &blanketVertices,
            const EdgeSet &blanketEdges);

    TopologyProvider *chooseTopologyProvider();

    void computeTargetInformation(
            const VertexSet &toRemove);

    JacobianMapping buildJacobianMapping(
            const g2o::OptimizableGraph::EdgeContainer &edges);

    void updateInputGraph(
            const VertexSet &toRemove,
            const EdgeSet &blanketEdges,
            const g2o::OptimizableGraph::EdgeContainer &newedges,
            const std::list<g2o::MatrixXD> &infomats);

    void cleanup();

    VertexSet mapForward(const VertexSet &vset);


private:
    SparsityOptions _opts;
    g2o::OptimizableGraph *_graph;
    EdgeMap *_edgeLookup;

    g2o::SparseOptimizer *_subgraph;            // d
    VertexMapping _forwardMap, _backwardMap;       // c
    g2o::MatrixXD _information;                 // c
    g2o::OptimizableGraph::VertexContainer _indexMapping;   // c
    std::list<TopologyProvider *> _topologies;          // d
    TopologyProvider::VertexInfoSet _vertexTypes;       // c
    TopologyProvider::EdgeInfoSet _edgeTypes;           // c
};

#endif /* VERTEX_REMOVER_H_ */
