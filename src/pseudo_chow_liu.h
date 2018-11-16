/*
 * pseudo_chow_liu.h
 *
 *  Created on: Dec 18, 2013
 *      Author: mazuran
 */

#ifndef PSEUDO_CHOW_LIU_H_
#define PSEUDO_CHOW_LIU_H_

#include <Eigen/Core>
#include <g2o/core/optimizable_graph.h>
#include <queue>
#include <list>
#include "sparsity_options.h"

class PseudoChowLiu {
public:
    typedef std::list<std::pair<g2o::OptimizableGraph::Vertex *,
            g2o::OptimizableGraph::Vertex *>> CorrelatedSkeletonTree;
    typedef std::list<CorrelatedSkeletonTree> SparsityPattern;

    PseudoChowLiu();
    virtual ~PseudoChowLiu();

    void setSparsityOptions(const SparsityOptions &opts) { _opts = opts; }
    void setTargetInformation(const Eigen::MatrixXd &info) { _information = &info; }
    void setVertices(const g2o::OptimizableGraph::VertexContainer &vertices) { _vertices = vertices; }

    void computeSparsityPattern();
    const SparsityPattern &getSparsityPattern() const { return _pattern; }

    Eigen::MatrixXd marginal(g2o::OptimizableGraph::Vertex *vert) const;
    Eigen::MatrixXd marginal(int vert) const;
    Eigen::MatrixXd jointMarginal(
            g2o::OptimizableGraph::Vertex *vert1,
            g2o::OptimizableGraph::Vertex *vert2) const;
    Eigen::MatrixXd jointMarginal(int vert1, int vert2) const;

private:
    void computeLookup();
    double weight(int vert1, int vert2);
    void fillEdges();
    void doKruskal();
    void fillCliques();
    Eigen::MatrixXd marginal(const std::vector<int> &keep) const;
    int vertexId(g2o::OptimizableGraph::Vertex *v) const;

    struct WeightedEdge {
        double weight;
        int vert1, vert2;
        bool operator<(const WeightedEdge &e) const { return weight < e.weight; }
    };

private:
    SparsityOptions _opts;
    const Eigen::MatrixXd *_information;
    Eigen::MatrixXd _pseudoCovariance;
    std::priority_queue<WeightedEdge> _edges;
    std::vector<WeightedEdge> _edgeBin;
    g2o::OptimizableGraph::VertexContainer _vertices;
    std::vector<int> _lookup;
    SparsityPattern _pattern;
};

#endif /* PSEUDO_CHOW_LIU_RETRIEVER_H_ */
