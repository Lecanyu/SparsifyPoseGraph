/*
 * vertex_remover.cpp
 *
 *  Created on: Oct 16, 2014
 *      Author: mazuran
 */

#include "vertex_remover.h"
#include "info_block_solver.h"
#include "type_info.h"
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include "optimizer.h"
#include "utils.h"

#if 0
std::ostream &operator<<(std::ostream &s, const g2o::HyperGraph::Vertex *v)
{
    return s << (void *) v << "[" << v->id() << "]";
}

std::ostream &operator<<(std::ostream &s, const g2o::HyperGraph::Edge *e)
{
    return s << e->vertices();
}
#endif

typedef InfoBlockSolver<g2o::BlockSolverTraits<-1, -1> > BlockSolver;
typedef g2o::LinearSolverCholmod<BlockSolver::PoseMatrixType> LinearSolver;
typedef g2o::OptimizationAlgorithmLevenberg OptimizationAlgorithm;


VertexRemover::VertexRemover() : _graph(NULL), _subgraph(new g2o::SparseOptimizer)
{
    LinearSolver *linearSolver = new LinearSolver;
    linearSolver->setBlockOrdering(false);
    BlockSolver *blockSolver = new BlockSolver(linearSolver);
    OptimizationAlgorithm *solver = new OptimizationAlgorithm(blockSolver);
    _subgraph->setAlgorithm(solver);
}

VertexRemover::~VertexRemover()
{
    delete _subgraph;
    for(TopologyProvider *topology: _topologies) {
        delete topology;
    }
}

void VertexRemover::setSparsityOptions(const SparsityOptions &opts)
{
    _opts = opts;
    for(std::list<TopologyProvider *>::const_iterator it = _topologies.begin();
            it != _topologies.end(); ++it) {
        (*it)->setSparsityOptions(opts);
    }
}

static bool vertexCompare(g2o::OptimizableGraph::Vertex *a, g2o::OptimizableGraph::Vertex *b) {
    return a->id() < b->id();
}

static bool pointerCompare(g2o::OptimizableGraph::Vertex *a, g2o::OptimizableGraph::Vertex *b) {
    return a < b;
}

VertexRemover::VertexSet VertexRemover::mapForward(const VertexSet &vset)
{
    VertexSet ret(vertexCompare);
    for(auto v: vset) {
        ret.insert(_forwardMap[v]);
    }
    return ret;
}

g2o::OptimizableGraph::EdgeContainer VertexRemover::remove(
        g2o::OptimizableGraph::Vertex *toRemove) {
    std::vector<g2o::OptimizableGraph::Vertex *> rmVertices;
    rmVertices.push_back(toRemove);
    return remove(rmVertices);
}

g2o::OptimizableGraph::EdgeContainer VertexRemover::remove(
        const std::vector<g2o::OptimizableGraph::Vertex *> &toRemove)
{
    g2o::OptimizableGraph::EdgeContainer addedEdges;VertexSet toRemoveSet(pointerCompare), deleted(pointerCompare);
    toRemoveSet.insert(toRemove.begin(), toRemove.end());

    for(size_t i = 0; i < toRemove.size(); i++) {
        VertexSet vmarkov, toRemoveNow(vertexCompare);
        if(deleted.count(toRemove[i]) > 0) continue;

        if(_opts.topology == SparsityOptions::Dense || _opts.topology == SparsityOptions::CliqueyDense) {
            vmarkov = extendedMarkovBlanketVertices(toRemove[i], toRemoveSet, toRemoveNow);
        } else {
            vmarkov = markovBlanketVertices(toRemove[i]);
            toRemoveNow.insert(toRemove[i]);
        }

        EdgeSet emarkov = markovBlanketEdges(vmarkov, toRemoveNow);
//        for(auto it:emarkov)
//        {
//            std::cout<<"edge "<<(*it).id()<<", ";
//            for(auto v:(*it).vertices())
//                std::cout<<(*v).id()<<"\n";
//        }

        buildSubgraph(toRemoveNow, vmarkov, emarkov);
        computeTargetInformation(mapForward(toRemoveNow));
        TopologyProvider *tp = chooseTopologyProvider();
        g2o::OptimizableGraph::EdgeContainer edges = tp->topology(_information, _indexMapping);
        std::list<g2o::MatrixXD> infos;

//        for(auto it:edges)
//        {
//            std::cout<<"edge "<<(*it).id()<<", ";
//            for(auto v:(*it).vertices())
//                std::cout<<(*v).id()<<" ";
//            std::cout<<"\n";
//        }

        if(edges.size() > 0) {
            if(tp->requiresOptimization()) {
                JacobianMapping mapping = buildJacobianMapping(edges);
                infos = optimizeInformation(mapping, _information);
            } else {
                for(g2o::OptimizableGraph::Edge *e: edges) {
                    Eigen::Map<g2o::MatrixXD> info(e->informationData(), e->dimension(), e->dimension());
                    infos.push_back(info);
                }
            }
        }

        updateInputGraph(toRemoveNow, emarkov, edges, infos);
        cleanup();
        addedEdges.insert(addedEdges.end(), edges.begin(), edges.end());
        deleted.insert(toRemoveNow.begin(), toRemoveNow.end());
    }
    return addedEdges;
}

VertexRemover::VertexSet VertexRemover::extendedMarkovBlanketVertices(
        g2o::OptimizableGraph::Vertex *root, const VertexSet &pickBin,
        VertexSet &picked) const
{
    picked = VertexSet(vertexCompare);
    VertexSet ret = markovBlanketVertices(root);

    picked.insert(root);
#if 0
    std::map<g2o::OptimizableGraph::Vertex *, VertexSet> allBlankets;

    for(auto v: pickBin) {
        if(v != root) {
            allBlankets[v] = markovBlanketVertices(v);
        }
    }
    bool added = true;
    while(added) {
        added = false;
        std::map<g2o::OptimizableGraph::Vertex *, VertexSet>::iterator it = allBlankets.begin();
        while(it != allBlankets.end()) {
            if(_opts.includeIntraClique) {
                if(setIntersection(ret, it->second).size() > 0) {
                    picked.insert(it->first);
                    ret.insert(it->second.begin(), it->second.end());
                    allBlankets.erase(it++);
                    added = true;
                } else {
                    ++it;
                }
            } else {
                if(ret.count(it->first) > 0) {
                    picked.insert(it->first);
                    ret.insert(it->second.begin(), it->second.end());
                    allBlankets.erase(it++);
                    added = true;
                } else {
                    ++it;
                }
            }
        }
    }
#else
    for(auto v: ret) {
        if(pickBin.count(v) > 0 && picked.count(v) == 0) {
            picked.insert(v);
            VertexSet other = markovBlanketVertices(v);
            ret.insert(other.begin(), other.end());
        }
    }
#endif

    return ret;
}

VertexRemover::VertexSet VertexRemover::markovBlanketVertices(
        g2o::OptimizableGraph::Vertex *root) const
{
    VertexSet vset(vertexCompare);

    vset.insert(root);

    /* Lookup all nodes in Markov blanket */
    for(EdgeSet::const_iterator ite = root->edges().begin();
            ite != root->edges().end(); ++ite) {
        const g2o::HyperGraph::VertexContainer &vertices = (*ite)->vertices();
        for(g2o::HyperGraph::VertexContainer::const_iterator itv = vertices.begin();
                itv != vertices.end(); ++itv) {
            vset.insert(static_cast<g2o::OptimizableGraph::Vertex *>(*itv));
        }
    }

    return vset;
}

VertexRemover::EdgeSet VertexRemover::markovBlanketEdges(
        g2o::OptimizableGraph::Vertex *root) const
{
    VertexSet hub(vertexCompare);
    hub.insert(root);
    return markovBlanketEdges(markovBlanketVertices(root), hub);
}

VertexRemover::EdgeSet VertexRemover::markovBlanketEdges(
        const VertexSet &mbVertices, const VertexSet &hubs) const
{
    EdgeSet edges;
    for(auto &vertex: mbVertices) {
        /* For each edge, if all endpoints are in the Markov blanket add it to internalEdges */
        for(auto &edge: vertex->edges()) {
            bool is_markov = true, found_hub = false;
            for(auto &vert: edge->vertices()) {
                g2o::OptimizableGraph::Vertex *v = static_cast<
                        g2o::OptimizableGraph::Vertex *>(vert);
                if(mbVertices.count(v) == 0) {
                    is_markov = false;
                    break;
                }

                if(hubs.count(v) > 0) {
                    found_hub = true;
                }
            }
            if(is_markov && (_opts.includeIntraClique || found_hub)) {
                edges.insert(static_cast<g2o::OptimizableGraph::Edge *>(edge));
            }
        }
    }
    return edges;
}

void VertexRemover::addSubgraphVertex(
        int id, g2o::OptimizableGraph::Vertex *original,
        const g2o::OptimizableGraph::Vertex *reparam)
{
    const VertexInfo *vi = typeinfo(original);
    assert(vi != typeinfo.vertexUnknown() && "Unknown type of vertex");
    g2o::OptimizableGraph::Vertex *newv = static_cast<
            g2o::OptimizableGraph::Vertex *>(vi->clone(original));
    _vertexTypes.insert(vi);
    newv->setId(id);
    _forwardMap[original] = newv;
    _backwardMap[newv] = original;
    if(reparam)
        vi->reparametrize(reparam, newv);
    _subgraph->addVertex(newv);
}

void VertexRemover::addSubgraphEdge(const g2o::HyperGraph::Edge *e)
{
    const EdgeInfo *ei = typeinfo(e);
    g2o::HyperGraph::Edge *clone = ei->clone(e);
    _edgeTypes.insert(ei);
    assert(ei != typeinfo.edgeUnknown() && "Unknown type of edge");
    for(auto &v: clone->vertices()) {
        g2o::OptimizableGraph::Vertex *vv = static_cast<
                g2o::OptimizableGraph::Vertex *>(v);
        assert(_forwardMap.count(vv) > 0 && "Unknown vertices in cloned edge");
        v = _forwardMap[vv];
    }
    _subgraph->addEdge(clone);
}

void VertexRemover::buildSubgraph(
            const VertexSet &toRemove,
            const VertexSet &blanketVertices,
            const EdgeSet &blanketEdges)
{
    int k = 0;
    assert(toRemove.size() > 0 && blanketVertices.size() > 0 && blanketEdges.size() > 0);

    /*
     * TODO: temporarily disabled reparametrization, looks like it gives slightly worse
     * kld values, although I'm not too sure why. The information matrix should be merely
     * rotated. Dunno.
     */
    /* Check if we can reparamerize all the vertices wrt toRemove */
    bool canReparametrize = false, closedFormEstimate = false;
    for(auto &v: blanketVertices) {
        canReparametrize = canReparametrize && typeinfo(v)->canReparametrize(*toRemove.begin());
    }

    if(_opts.linPoint == SparsityOptions::Local) {

        /*
         * Local linearization point. In order to compute a closed form solution
         * the following needs to hold:
         *  - all vertices (except toRemove) appear only once in the edges
         *  - initialEstimatePossible() works for all edges with toRemove fixed
         */

        std::map<g2o::HyperGraph::Vertex *, int> nconnections;
        std::set<g2o::HyperGraph::Vertex *> fixed;
        fixed.insert(*toRemove.begin());

        closedFormEstimate = true;
        for(auto v: blanketVertices) {
            if(v != *toRemove.begin())
                nconnections[v] = 0;
        }

        for(auto e: blanketEdges) {
            g2o::OptimizableGraph::Edge *ee = static_cast<
                        g2o::OptimizableGraph::Edge *>(e);
            for(auto v: e->vertices()) {
                g2o::OptimizableGraph::Vertex *vv = static_cast<
                            g2o::OptimizableGraph::Vertex *>(v);
                if(v != *toRemove.begin()) {
                    nconnections[v]++;
                    closedFormEstimate = closedFormEstimate &&
                            ee->initialEstimatePossible(fixed, vv) > 0;
                }
            }
        }

        for(auto &nconn: nconnections) {
            if(nconn.second > 1) {
                closedFormEstimate = false;
                break;
            }
        }
    }

    /*
     * Add all the vertices, and if possible reparametrize them all so that
     * the one to be removed is the origin.
     */
    for(auto v: toRemove) {
        addSubgraphVertex(k++, v, canReparametrize ? *toRemove.begin() : NULL);
    }

    for(auto v: blanketVertices) {
        if(toRemove.count(v) == 0)
            addSubgraphVertex(k++, v, canReparametrize ? *toRemove.begin() : NULL);
    }

    /* Copy all of the edges */
    for(auto e: blanketEdges) {
        addSubgraphEdge(static_cast<const g2o::OptimizableGraph::Edge *>(e));
    }

    if(closedFormEstimate) {
        /*
         * Local linearization point. Compute the closed form solution with
         * initialEstimate()
         */
        std::set<g2o::HyperGraph::Vertex *> fixed;
        fixed.insert(_forwardMap[*toRemove.begin()]);
        typeinfo(*toRemove.begin())->setZero(_forwardMap[*toRemove.begin()]);
        for(auto e: _subgraph->edges()) {
            g2o::OptimizableGraph::Edge *ee = static_cast<
                        g2o::OptimizableGraph::Edge *>(e);
            for(auto v: e->vertices()) {
                g2o::OptimizableGraph::Vertex *vv = static_cast<
                            g2o::OptimizableGraph::Vertex *>(v);
                if(v != *toRemove.begin()) {
                    ee->initialEstimate(fixed, vv);
                }
            }
        }
    } else if(_opts.linPoint != SparsityOptions::Global) {
        /*
         * Local linearization point without closed form.
         * Optimize the graph by fixing the node to remove.
         */
        _forwardMap[*toRemove.begin()]->setFixed(true);
        _subgraph->initializeOptimization();
        _subgraph->optimize(10);
        _forwardMap[*toRemove.begin()]->setFixed(false);
    }
}

void VertexRemover::computeTargetInformation(const VertexSet &toRemove)
{
    assert(toRemove.size() > 0);
    _subgraph->initializeOptimization();
    OptimizationAlgorithm *solver = static_cast<OptimizationAlgorithm *>(_subgraph->solver());
    BlockSolver *blockSolver = static_cast<BlockSolver *>(solver->solver());
    solver->init(false);
    blockSolver->buildStructure();
    blockSolver->buildSystem();

    int removalSize = 0;
    for(auto v: toRemove) {
        removalSize += v->dimension();
    }

    Eigen::SparseMatrix<double> info = sparse(blockSolver->information());
    std::vector<int> indicesKeep, indicesMarginalize, verticesMarginalize;

    indicesKeep.reserve(info.cols() - removalSize);
    indicesMarginalize.reserve(removalSize);
    verticesMarginalize.reserve(toRemove.size());
    _indexMapping.reserve(_subgraph->vertices().size() - toRemove.size());

    int margidx = 0, k = 0, i = 0;

    for(auto *v: _subgraph->indexMapping()) {
        if(toRemove.count(v) > 0) {
            verticesMarginalize.push_back(i);
            for(int j = 0; j < v->dimension(); j++) {
                indicesMarginalize.push_back(k + j);
            }
        } else {
            for(int j = 0; j < v->dimension(); j++) {
                indicesKeep.push_back(k + j);
            }
        }
        k += v->dimension();
        i++;
    }

    std::sort(verticesMarginalize.begin(), verticesMarginalize.end());
    for(int i = 0, j = 0; i < _subgraph->indexMapping().size(); i++) {
        if(j < verticesMarginalize.size() && verticesMarginalize[j] == i) {
            j++;
        } else {
            _indexMapping.push_back(_subgraph->indexMapping()[i]);
        }
    }

    /* Compute Schur complement */
    Eigen::LLT<g2o::MatrixXD> chol(selectVariables(info, indicesMarginalize));
    Eigen::SparseMatrix<double> mixed = selectVariables(info, indicesMarginalize, indicesKeep);
    _information = selectVariables(info, indicesKeep);
    _information -= mixed.transpose() * chol.solve(g2o::MatrixXD(mixed));
    _information.triangularView<Eigen::StrictlyLower>() =
            _information.triangularView<Eigen::StrictlyUpper>().transpose();
}

TopologyProvider *VertexRemover::chooseTopologyProvider()
{
    TopologyProvider *chosen = NULL;
    for(TopologyProvider *tp: _topologies) {
        if(tp->applicable(_vertexTypes, _edgeTypes)) {
            chosen = tp;
            break;
        }
    }
    assert(chosen && "No valid topology provider for Markov blanket");
    return chosen;
}


JacobianMapping VertexRemover::buildJacobianMapping(
        const g2o::OptimizableGraph::EdgeContainer &edges)
{
    std::map<g2o::OptimizableGraph::Vertex *, int> lookup;
    JacobianMapping mapping;
    g2o::JacobianWorkspace jw;
    int k = 0;

    for(g2o::OptimizableGraph::Vertex *v: _indexMapping) {
        lookup[v] = k;
        k += v->dimension();
    }

    for(g2o::OptimizableGraph::Edge *e: edges) {
        jw.updateSize(e);
    }

    jw.allocate();

    for(g2o::OptimizableGraph::Edge *e: edges) {
        e->linearizeOplus(jw);
        mapping.push_back(MeasurementJacobian());
        for(int i = 0; i < e->vertices().size(); i++) {
            g2o::OptimizableGraph::Vertex *v = static_cast<
                    g2o::OptimizableGraph::Vertex *>(e->vertex(i));
            Eigen::Map<g2o::MatrixXD> J(jw.workspaceForVertex(i),
                    e->dimension(), v->dimension());
            mapping.back().push_back(std::make_pair(g2o::MatrixXD(J), lookup[v]));
        }
    }

    return mapping;
}

void VertexRemover::updateInputGraph(
        const VertexSet &toRemove,
        const EdgeSet &blanketEdges,
        const g2o::OptimizableGraph::EdgeContainer &newedges,
        const std::list<g2o::MatrixXD> &infomats)
{
//    std::cout<<"remove edge: \n";
    for(g2o::HyperGraph::Edge *e: blanketEdges) {
//        for(auto v:e->vertices())
//            std::cout<<(*v).id()<<", ";
//        std::cout<<"\n";
        _graph->removeEdge(e);

        /**
         * Canyu Le
         * To fix memory leak, we delete the _edgeLookup here to guarantee consistent result.
         * */
        // -----------------------
        _edgeLookup->erase(e);
        // -----------------------
    }

//    std::cout<<"remove vertices: \n";
    for(g2o::HyperGraph::Vertex *v: toRemove) {
//        std::cout<<v->id()<<"\n";
        _graph->removeVertex(v);
    }

//    std::cout<<"add new edges: \n";
    std::list<g2o::MatrixXD>::const_iterator it = infomats.begin();
    for(g2o::OptimizableGraph::Edge *e: newedges) {
        Eigen::Map<g2o::MatrixXD> info(
                e->informationData(), e->dimension(), e->dimension());
        info = *it; ++it;
        for(int i = 0; i < e->vertices().size(); i++) {
            g2o::OptimizableGraph::Vertex *v = static_cast<
                    g2o::OptimizableGraph::Vertex *>(e->vertex(i));
            e->setVertex(i, _backwardMap[v]);
        }
        _graph->addEdge(e);

//        for(auto v:e->vertices())
//            std::cout<<(*v).id()<<", ";
//        std::cout<<"\n";
    }
    std::cout<<"";
}

void VertexRemover::cleanup()
{
    _subgraph->clear();
    _forwardMap.clear();
    _backwardMap.clear();
    _vertexTypes.clear();
    _edgeTypes.clear();
    _indexMapping.clear();
    _information.resize(0, 0);
}

