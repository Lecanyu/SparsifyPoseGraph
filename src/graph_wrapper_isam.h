/*
 * graph_wrapper_isam.h
 *
 *  Created on: 15/dic/2013
 *      Author: Mladen Mazuran
 */

#ifndef GRAPH_WRAPPER_ISAM_H_
#define GRAPH_WRAPPER_ISAM_H_

#include "graph_wrapper.h"
#include <isam/Slam.h>
#include <map>
#include <Eigen/Sparse>
#include "auto_delete_map.h"

class GraphWrapperISAM : public GraphWrapper
{
public:
    class Vertex;
    class Edge;

    typedef AutoDeleteMap<const isam::Node *,   GraphWrapper::Vertex *> VertexMap;
    typedef AutoDeleteMap<const isam::Factor *, GraphWrapper::Edge *>   EdgeMap;

    class Vertex : public GraphWrapper::Vertex {
    public:
        Vertex(int id, isam::Node *v, const EdgeMap &l);
        std::vector<const GraphWrapper::Edge *> edges() const;
        int id() const { return _id; }
        IsometryXd estimate() const;
        bool is2d() const;
    private:
        int _id;
        isam::Node *_v;
        const EdgeMap *_l;
    };

    class Edge : public GraphWrapper::Edge {
    public:
        Edge(isam::Factor *e, const VertexMap &l);
        std::vector<const GraphWrapper::Vertex *> vertices() const;
        IsometryXd measurement() const;
        Eigen::MatrixXd information() const;
    private:
        isam::Factor *_e;
        const VertexMap *_l;
    };

public:
    GraphWrapperISAM();
    virtual ~GraphWrapperISAM();

    void addVertex(int id, const IsometryXd &init);
    void addEdge(int from, int to, const IsometryXd &meas, const Eigen::MatrixXd &info);
    void addEdge(int to, const IsometryXd &meas, const Eigen::MatrixXd &info);

    void addEdgeSqrt(int from, int to, const IsometryXd &meas, const Eigen::MatrixXd &sqrtinfo);
    void addEdgeSqrt(int to, const IsometryXd &meas, const Eigen::MatrixXd &sqrtinfo);

    void optimize();

    GraphWrapperISAM *clonePortion(int maxid);

    Eigen::VectorXd estimate();
    Eigen::MatrixXd information();
    Eigen::MatrixXd covariance();
    
    void marginalize(const std::vector<int> &which, const SparsityOptions &options);
    
    GraphWrapper::Vertex *vertex(int id) { return _vertexLookup[_forwardLookup[id]]; }

    void write(std::ofstream &s);

    double kullbackLeibler(GraphWrapper *other);

    double chi2(GraphWrapper *other) { return nan(""); }

    double chi2() const { return _slam->chi2(); }

    void debugPrint(std::ostream &s) const;

    void printStats(std::ostream &s) const;

    Eigen::VectorXd estimateDifference(GraphWrapper *other);

    void setEstimate(int vertexid, const IsometryXd &est);

    void saveInformation(const char *fname);

private:
    isam::Slam *_slam;
    std::map<int, isam::Node *> _forwardLookup;
    std::map<isam::Node *, int> _reverseLookup;
    VertexMap _vertexLookup;
    EdgeMap _edgeLookup;
};

#endif /* GRAPH_WRAPPER_ISAM_H_ */
