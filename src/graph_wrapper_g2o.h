/*
 * graph_wrapper_g2o.h
 *
 *  Created on: 15/dic/2013
 *      Author: Mladen Mazuran
 */

#ifndef GRAPH_WRAPPER_G2O_H_
#define GRAPH_WRAPPER_G2O_H_

#include "graph_wrapper.h"
#include <g2o/core/sparse_optimizer.h>
#include <Eigen/Sparse>
#include "se2_compatibility.h"
#include "se3_compatibility.h"
#include "auto_delete_map.h"

class GraphWrapperISAM;

typedef VertexSE2ISAM VertexSE2Type;
typedef VertexSE3ISAM VertexSE3Type;
//typedef g2o::VertexSE3 VertexSE3Type;

typedef EdgeSE2ISAM   EdgeSE2Type;
typedef EdgeSE3ISAM   EdgeSE3Type;
//typedef g2o::EdgeSE3 EdgeSE3Type;

class GraphWrapperG2O : public GraphWrapper
{
public:
    class Vertex;
    class Edge;

    typedef AutoDeleteMap<const g2o::HyperGraph::Vertex *, GraphWrapper::Vertex *> VertexMap;
    typedef AutoDeleteMap<const g2o::HyperGraph::Edge *,   GraphWrapper::Edge *>   EdgeMap;

    class Vertex : public GraphWrapper::Vertex {
    public:
        Vertex(g2o::HyperGraph::Vertex *v, const EdgeMap &l);
        std::vector<const GraphWrapper::Edge *> edges() const;
        int id() const { return _v->id(); }
        IsometryXd estimate() const;
        bool is2d() const;
    private:
        g2o::HyperGraph::Vertex *_v;
        const EdgeMap *_l;
    };

    class Edge : public GraphWrapper::Edge {
    public:
        Edge(g2o::HyperGraph::Edge *e, const VertexMap &l);
        std::vector<const GraphWrapper::Vertex *> vertices() const;
        IsometryXd measurement() const;
        Eigen::MatrixXd information() const;
    private:
        g2o::HyperGraph::Edge *_e;
        const VertexMap *_l;
    };

public:        
    GraphWrapperG2O(bool verbose = false, bool useGLC = false);
    GraphWrapperG2O(const char *fname, bool optimizeit = true, bool useGLC = false);
    virtual ~GraphWrapperG2O();

    GraphWrapperISAM *toISAM();

    void addVertex(int id, const IsometryXd &init);
    void addEdge(int from, int to, const IsometryXd &meas, const Eigen::MatrixXd &info);

    void optimize();
    void optimizeWithReference(GraphWrapperG2O *other);

    GraphWrapperG2O *clonePortion(int maxid);

    Eigen::VectorXd estimate();
    Eigen::MatrixXd information();
    Eigen::MatrixXd covariance();

    Eigen::SparseMatrix<double> sparseInformation() const;

    GraphWrapper::Vertex *vertex(int id);

    void marginalize(const std::vector<int> &which, const SparsityOptions &options);
    void marginalizeNoOptimize(const std::vector<int> &which, const SparsityOptions &options);

    void write(std::ofstream &s);

    double kullbackLeibler(GraphWrapper *other);

    double chi2(GraphWrapper *other);

    double chi2() const { _so->computeActiveErrors(); return _so->chi2(); }

    void debugPrint(std::ostream &s) const;

    void saveInformation(const char *fname);

    void printStats(std::ostream &s) const;

    Eigen::VectorXd estimateDifference(GraphWrapper *other);

    void setEstimate(int vertexid, const IsometryXd &est);

    void setRobust(bool robust) { _robust = robust; }

    void push();
    void pop();

    /**
     * Canyu Le
     * Fix a specific vertex to optimize instead of default 0.
     * */
    // -----------------------
    void optimizeFromId(int first_vertex_id);
    // -----------------------

    /**
     * Canyu Le
     * For debug and checking graph status.
     * */
    // -----------------------
    std::string getStats()
    {
        std::stringstream ss;
        Eigen::SparseMatrix<double> info = sparseInformation();
        double fillin = info.nonZeros() / double(info.rows() * info.cols());
        ss << "nodes = " << _so->vertices().size() - 1 << "; edges = " << _so->edges().size() <<
          "; fillin = " << fillin * 100 << "%";
        return ss.str();
    }
    // -----------------------

    /**
     * Canyu Le
     * Write some get methods to access private variable.
     * */
    // -----------------------
    g2o::SparseOptimizer* getOptimizer(){ return _so; }
    VertexMap getVertexMap(){ return _vertexLookup; }
    EdgeMap getEdgeMap(){ return _edgeLookup; }
    // -----------------------

private:
    g2o::SparseOptimizer *newSparseOptimizer(bool verbose = false) const;

    void computeIndices(
            GraphWrapper *other,
            std::vector<int> &indicesKeep,
            std::vector<int> &indicesMarginalize);

private:
    g2o::SparseOptimizer *_so;
    VertexMap _vertexLookup;
    EdgeMap _edgeLookup;
    bool _robust, _verbose, _useGLC;
};

#endif /* GRAPH_WRAPPER_G2O_H_ */
