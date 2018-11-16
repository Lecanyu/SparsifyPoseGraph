/*
 * pseudo_chow_liu.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: mazuran
 */

#include "pseudo_chow_liu.h"
#include <Eigen/Cholesky>
#include <set>
#include "utils.h"

PseudoChowLiu::PseudoChowLiu() : _information(NULL)
{
}

PseudoChowLiu::~PseudoChowLiu()
{
}

void PseudoChowLiu::computeLookup()
{
    int k = 0;
    _lookup.clear();
    _lookup.push_back(0);
    for(g2o::OptimizableGraph::VertexContainer::const_iterator it = _vertices.begin();
            it != _vertices.end(); ++it) {
        k += (*it)->dimension();
        _lookup.push_back(k);
    }
}

void PseudoChowLiu::computeSparsityPattern()
{
    computeLookup();

    assert(_information && "Target information matrix not set");
    assert(_information->cols() == _lookup.back() &&
            "Discrepancy between information matrix and vertices");

    _pattern.clear();
    int n = _vertices.size();
    int m = int((1 + _opts.chordRatio) * (n - 1));
    bool full = (m >= n * (n - 1) / 2);

    if(n == 2) {
        /* No need to compute anything for only two nodes */
        _pattern.push_back(CorrelatedSkeletonTree(1,
                std::make_pair(_vertices[0], _vertices[1])));
    } else if(_opts.topology == SparsityOptions::Dense || (
            _opts.topology == SparsityOptions::Subgraph && full)) {
        /* Dense graph with uncorrelated edges */
        for(int i = 0; i < n - 1; i++) {
            for(int j = i + 1; j < n; j++) {
                _pattern.push_back(CorrelatedSkeletonTree(1,
                        std::make_pair(_vertices[i], _vertices[j])));
            }
        }
    } else {

        fillEdges();
        doKruskal();

        if(_opts.topology == SparsityOptions::Tree ||
                _opts.topology == SparsityOptions::Subgraph) {
            /* Tree/subgraph of uncorrelated edges */
            int nedges = (_opts.topology == SparsityOptions::Tree) ? n - 1 : m;

            for(int i = 0; i < nedges; i++) {
                _pattern.push_back(CorrelatedSkeletonTree(1, std::make_pair(
                        _vertices[_edgeBin[i].vert1], _vertices[_edgeBin[i].vert2])));
            }
        } else if(_opts.topology == SparsityOptions::CliqueyDense || (
                _opts.topology == SparsityOptions::CliqueySubgraph && full)) {
            /* Tree of fully correlated edges */
            CorrelatedSkeletonTree tree;
            for(int i = 0; i < n - 1; i++) {
                tree.push_back(std::make_pair(
                        _vertices[_edgeBin[i].vert1], _vertices[_edgeBin[i].vert2]));
            }
            _pattern.push_back(tree);
        } else /* if(_opts.topology == SparsityOptions::CliqueySubgraph) */ {
            /* Tree of partially correlated edges */
            fillCliques();
        }
    }
}

static std::vector<int> range(int min, int lessthan)
{
    std::vector<int> ret;
    for(; min < lessthan; min++) {
        ret.push_back(min);
    }
    return ret;
}

static std::vector<int> join(const std::vector<int> &a, const std::vector<int> &b)
{
    std::vector<int> ret;
    ret.insert(ret.end(), a.begin(), a.end());
    ret.insert(ret.end(), b.begin(), b.end());
    return ret;
}

static std::vector<int> complement(const std::vector<int> &a, int bound)
{
    std::vector<int> ret;
    for(int i = 0, j = 0; i < bound; i++) {
        if(j < int(a.size()) && a[j] == i) {
            j++;
        } else {
            ret.push_back(i);
        }
    }
    return ret;
}

int PseudoChowLiu::vertexId(g2o::OptimizableGraph::Vertex *v) const
{
    for(size_t i = 0; i < _vertices.size(); i++) {
        if(_vertices[i] == v) {
            return i;
        }
    }
    assert(false && "vertex does not exist");
    return -1;
}

Eigen::MatrixXd PseudoChowLiu::marginal(const std::vector<int> &keep) const
{
    std::vector<int> marginalize = complement(keep, _information->rows());
    Eigen::MatrixXd mixed = selectVariables(*_information, keep, marginalize);
    Eigen::LLT<Eigen::MatrixXd> chol(selectVariables(*_information, marginalize));
    Eigen::MatrixXd schur = selectVariables(*_information, keep) -
            mixed * chol.solve(mixed.transpose());
    return schur.selfadjointView<Eigen::Upper>();
}

Eigen::MatrixXd PseudoChowLiu::marginal(g2o::OptimizableGraph::Vertex *vert) const
{
    return marginal(vertexId(vert));
}

Eigen::MatrixXd PseudoChowLiu::marginal(int vert) const
{
    const int d = static_cast<const g2o::OptimizableGraph::Vertex *>(_vertices[vert])->dimension();
    std::vector<int> keep = range(_lookup[vert], _lookup[vert] + d);
    return marginal(keep);
}

Eigen::MatrixXd PseudoChowLiu::jointMarginal(
        g2o::OptimizableGraph::Vertex *vert1,
        g2o::OptimizableGraph::Vertex *vert2) const
{
    return jointMarginal(vertexId(vert1), vertexId(vert2));
}

Eigen::MatrixXd PseudoChowLiu::jointMarginal(int vert1, int vert2) const
{
    const int d1 = static_cast<const g2o::OptimizableGraph::Vertex *>(_vertices[vert1])->dimension(),
            d2 = static_cast<const g2o::OptimizableGraph::Vertex *>(_vertices[vert2])->dimension();
    std::vector<int> keep = join(
            range(_lookup[vert1], _lookup[vert1] + d1),
            range(_lookup[vert2], _lookup[vert2] + d2));
    return marginal(keep);
}

double PseudoChowLiu::weight(int vert1, int vert2)
{
    const int d1 = static_cast<const g2o::OptimizableGraph::Vertex *>(_vertices[vert1])->dimension(),
            d2 = static_cast<const g2o::OptimizableGraph::Vertex *>(_vertices[vert2])->dimension();

    Eigen::MatrixXd jointCov = selectVariables(_pseudoCovariance,
            join(range(_lookup[vert1], _lookup[vert1] + d1),
                    range(_lookup[vert2], _lookup[vert2] + d2)));

    Eigen::LDLT<Eigen::MatrixXd> x(_pseudoCovariance.block(_lookup[vert1], _lookup[vert1], d1, d1));
    Eigen::LDLT<Eigen::MatrixXd> y(_pseudoCovariance.block(_lookup[vert2], _lookup[vert2], d2, d2));
    Eigen::LDLT<Eigen::MatrixXd> xy(jointCov);

    return x.vectorD().array().log().sum() + y.vectorD().array().log().sum() - xy.vectorD().array().log().sum();
}

void PseudoChowLiu::fillEdges()
{
    static const double tikhonov_eps = 1;
    _edges = std::priority_queue<WeightedEdge>(); // Clear edges
    Eigen::LLT<Eigen::MatrixXd> llt(*_information + tikhonov_eps * Eigen::MatrixXd::Identity(_information->rows(), _information->cols()));
    _pseudoCovariance = llt.solve(Eigen::MatrixXd::Identity(_information->rows(), _information->cols()));
    for(int i = 0; i < (int) _vertices.size() - 1; i++) {
        for(int j = i + 1; j < (int) _vertices.size(); j++) {
            _edges.push({ weight(i, j), i, j });
        }
    }
}

void PseudoChowLiu::fillCliques()
{
    int n = _vertices.size();
    int m = int((1 + _opts.chordRatio) * (n - 1));

    std::vector<std::set<int> > cliques;
    cliques.reserve(n - 1);
    for(int i = 0; i < n - 1; i++) {
        cliques.push_back({ _edgeBin[i].vert1, _edgeBin[i].vert2 });
    }

    bool joined = true;
    int nedges = n - 1;

    for(int nedges = n - 1, maxfill = 1; nedges < m && joined; maxfill++) {
        joined = false;
        int minfill = std::numeric_limits<int>::max();
        for(int i = 0; i < cliques.size(); i++) {
            for(int j = i + 1; j < cliques.size(); j++) {
                if(setIntersection(cliques[i], cliques[j]).size() > 0) {
                    int thisfill = (cliques[i].size() - 1) * (cliques[j].size() - 1);
                    minfill = std::min(thisfill, minfill);
                    if(thisfill <= maxfill && nedges + thisfill <= m) {
                        cliques[i].insert(cliques[j].begin(), cliques[j].end());
                        nedges += thisfill;
                        cliques.erase(cliques.begin() + j);
                        joined = true;
                        j--;
                    }
                }
            }
        }
        if(!joined && minfill > maxfill) {
            joined = true;
            maxfill = minfill - 1;
        }
    }

    for(int i = 0; i < cliques.size(); i++) {
        _pattern.push_back(CorrelatedSkeletonTree());
    }

    for(int i = 0; i < n - 1; i++) {
        SparsityPattern::iterator it = _pattern.begin();
        for(int j = 0; j < cliques.size(); ++it, j++) {
            if(cliques[j].count(_edgeBin[i].vert1) > 0 &&
                    cliques[j].count(_edgeBin[i].vert2)) {
                it->push_back(std::make_pair(
                        _vertices[_edgeBin[i].vert1],
                        _vertices[_edgeBin[i].vert2]));
            }
        }
    }
}

void PseudoChowLiu::doKruskal()
{
    std::vector<std::set<int> > connectivity;
    std::vector<WeightedEdge> rejectBin, acceptBin;
    connectivity.reserve(_vertices.size());

    for(int i = 0; i < (int) _vertices.size(); i++) {
        connectivity.push_back(std::set<int>());
        connectivity.back().insert(i);
    }
    
    while(_edges.size() > 0) {
        int set1 = 0, set2 = 0;
        WeightedEdge e = _edges.top();
        _edges.pop();
        
        for(int i = 0; i < (int) connectivity.size(); i++) {
            if(connectivity[i].count(e.vert1) > 0) {
                set1 = i;
            }
            if(connectivity[i].count(e.vert2) > 0) {
                set2 = i;
            }
        }

        if(set1 != set2) {
            acceptBin.push_back(e);
            connectivity[set1].insert(connectivity[set2].begin(), connectivity[set2].end());
            connectivity.erase(connectivity.begin() + set2);
        } else {
            rejectBin.push_back(e);
        }
    }

    _edgeBin = acceptBin;
    _edgeBin.insert(_edgeBin.end(), rejectBin.begin(), rejectBin.end());
}

