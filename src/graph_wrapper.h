/*
 * graph_wrapper.h
 *
 *  Created on: 15/dic/2013
 *      Author: Mladen Mazuran
 */

#ifndef GRAPH_WRAPPER_H_
#define GRAPH_WRAPPER_H_

#include <Eigen/Core>
#include <vector>
#include <fstream>
#include "isometryxd.h"
#include "sparsity_options.h"

class GraphWrapper {
public: /* Types */
    class Edge;

    class Vertex {
    public:
        virtual ~Vertex() {}
        virtual int id() const = 0;
        virtual IsometryXd estimate() const = 0;
        virtual std::vector<const Edge *> edges() const = 0;
        virtual bool is2d() const = 0;
        virtual bool is3d() const { return !is2d(); }
    };

    class Edge {
    public:
        virtual ~Edge() {}
        virtual std::vector<const Vertex *> vertices() const = 0;
        virtual IsometryXd measurement() const = 0;
        virtual Eigen::MatrixXd information() const = 0;
    };

public: /* Methods */

    virtual ~GraphWrapper() {}

    virtual void addVertex(int id, const IsometryXd &init) = 0;

    virtual void addEdge(int from, int to, const IsometryXd &meas, const Eigen::MatrixXd &info) = 0;

    virtual void optimize() = 0;

    virtual GraphWrapper *clonePortion(int maxid) = 0;

    virtual Eigen::VectorXd estimate() = 0;

    virtual Eigen::MatrixXd information() = 0;

    virtual Eigen::MatrixXd covariance() = 0;

    virtual void marginalize(const std::vector<int> &which, const SparsityOptions &options) = 0;

    virtual void write(std::ofstream &s) = 0;

    virtual void write(const char *fname) { std::ofstream f(fname); write(f); f.close(); }

    virtual double kullbackLeibler(GraphWrapper *other) = 0;

    virtual double chi2(GraphWrapper *other) = 0;

    virtual double chi2() const = 0;

    virtual const std::vector<Vertex *> &vertices() const { return _vertices; }

    virtual Vertex *vertex(int id) = 0;

    virtual void debugPrint(std::ostream &s) const = 0;

    virtual void printStats(std::ostream &s) const = 0;

    virtual void setEstimate(int vertexid, const IsometryXd &est) = 0;

protected:
    std::vector<Vertex *> _vertices;
};

inline std::ostream &operator<<(std::ostream &s, const GraphWrapper &gw)
{
    gw.debugPrint(s);
    return s;
}

#endif /* GRAPH_WRAPPER_H_ */
