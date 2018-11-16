/*
 * graph_wrapper_g2o.cpp
 *
 *  Created on: 15/dic/2013
 *      Author: Mladen Mazuran
 */

#include "graph_wrapper_g2o.h"
#include <g2o/core/block_solver.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/robust_kernel_impl.h>
#include <g2o/types/slam2d/vertex_se2.h>
#include "info_block_solver.h"
#include "graph_wrapper_isam.h"
#include "utils.h"
#include "types.h"
#include "vertex_remover.h"
#include "topology_provider_se2.h"
#include "topology_provider_se3.h"
#include "topology_provider_glc.h"

typedef InfoBlockSolver<g2o::BlockSolverTraits<-1, -1> > SlamBlockSolver;
typedef g2o::LinearSolverCholmod<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;
typedef g2o::OptimizationAlgorithmLevenberg SlamOptimizationAlgorithm;

static bool vertexCompare(const GraphWrapper::Vertex *v, int id)
{
    return v->id() < id;
}

/* ---------------------------------------------------------------------------------------------- */
/*                                    GraphWrapperG2O::Vertex                                     */
/* ---------------------------------------------------------------------------------------------- */

GraphWrapperG2O::Vertex::Vertex(g2o::HyperGraph::Vertex *v, const GraphWrapperG2O::EdgeMap &l) :
        _v(v), _l(&l)
{
}

std::vector<const GraphWrapper::Edge *> GraphWrapperG2O::Vertex::edges() const
{
    std::vector<const GraphWrapper::Edge *> ret;
    for(const g2o::HyperGraph::Edge *e: _v->edges()) {
        ret.push_back(_l->at(e));
    }
    return ret;
}

IsometryXd GraphWrapperG2O::Vertex::estimate() const
{
    const VertexSE2Type *vse2 = dynamic_cast<const VertexSE2Type *>(_v);
    const VertexSE3Type *vse3 = dynamic_cast<const VertexSE3Type *>(_v);
    assert((vse2 || vse3) && "Can handle only SE2 or SE3 vertices for now");
    return vse2 ? IsometryXd(vse2->estimate()) : IsometryXd(vse3->estimate());
}

bool GraphWrapperG2O::Vertex::is2d() const
{
    return dynamic_cast<const VertexSE2Type *>(_v);
}

/* ---------------------------------------------------------------------------------------------- */
/*                                     GraphWrapperG2O::Edge                                      */
/* ---------------------------------------------------------------------------------------------- */

GraphWrapperG2O::Edge::Edge(g2o::HyperGraph::Edge *e, const GraphWrapperG2O::VertexMap &l) :
        _e(e), _l(&l)
{
}

std::vector<const GraphWrapper::Vertex *> GraphWrapperG2O::Edge::vertices() const
{
    std::vector<const GraphWrapper::Vertex *> ret;
    for(const g2o::HyperGraph::Vertex *v: _e->vertices()) {
        ret.push_back(_l->at(v));
    }
    return ret;
}

IsometryXd GraphWrapperG2O::Edge::measurement() const
{
    const EdgeSE2Type *ese2 = dynamic_cast<const EdgeSE2Type *>(_e);
    const EdgeSE3Type *ese3 = dynamic_cast<const EdgeSE3Type *>(_e);
    assert((ese2 || ese3) && "Can handle only SE2 or SE3 binary edges for now");
    return ese2 ? IsometryXd(ese2->measurement()) : IsometryXd(ese3->measurement());
}

Eigen::MatrixXd GraphWrapperG2O::Edge::information() const
{
    const EdgeSE2Type *ese2 = dynamic_cast<const EdgeSE2Type *>(_e);
    const EdgeSE3Type *ese3 = dynamic_cast<const EdgeSE3Type *>(_e);
    assert((ese2 || ese3) && "Can handle only SE2 or SE3 binary edges for now");
    return ese2 ? Eigen::MatrixXd(ese2->information()) : Eigen::MatrixXd(ese3->information());
}

/* ---------------------------------------------------------------------------------------------- */
/*                                        GraphWrapperG2O                                         */
/* ---------------------------------------------------------------------------------------------- */

GraphWrapperG2O::GraphWrapperG2O(bool verbose, bool useGLC) :
    _so(newSparseOptimizer(verbose)), _robust(false), _verbose(verbose), _useGLC(useGLC)
{
}

GraphWrapperG2O::GraphWrapperG2O(const char *fname, bool optimizeit, bool useGLC) :
    _so(newSparseOptimizer()), _robust(false), _verbose(false), _useGLC(useGLC)
{
    /* Load the file in a separate optimizer, so we can convert the nodes/edges to our format */
    g2o::SparseOptimizer *load = newSparseOptimizer(_verbose);
    load->load(fname);
    if(optimizeit) {
        load->vertex(0)->setFixed(true);
        load->initializeOptimization();
        load->optimize(50);
    }

    /* Create/convert nodes/edges */
    std::map<int, g2o::HyperGraph::Vertex *> vertices(load->vertices().begin(), load->vertices().end());
    for(auto v: vertices) {
        const g2o::VertexSE2 *vse2 = dynamic_cast<const g2o::VertexSE2 *>(v.second);
        const g2o::VertexSE3 *vse3 = dynamic_cast<const g2o::VertexSE3 *>(v.second);
        assert((vse2 || vse3) && "Can handle only SE2 or SE3 vertices for now");
        addVertex(v.first, vse2 ? IsometryXd(vse2->estimate()) : IsometryXd(vse3->estimate()));
    }

    for(auto e: load->edges()) {
        const g2o::EdgeSE2 *ese2 = dynamic_cast<const g2o::EdgeSE2 *>(e);
        const g2o::EdgeSE3 *ese3 = dynamic_cast<const g2o::EdgeSE3 *>(e);
        assert((ese2 || ese3) && "Can handle only SE2 or SE3 binary edges for now");
        IsometryXd meas;
        Eigen::MatrixXd info;
        if(ese3) {
            meas = IsometryXd(ese3->measurement());
#ifndef G2S_QUATERNIONS
            info = EdgeSE3ISAM::eulerInformationFromQuat(ese3->information());
#else /* G2S_QUATERNIONS */
            info = ese3->information();
#endif /* G2S_QUATERNIONS */
        } else {
            meas = IsometryXd(ese2->measurement());
            info = ese2->information();
        }
        addEdge(e->vertices()[0]->id(), e->vertices()[1]->id(), meas, info);
    }

    delete load;

    if(optimizeit) {
        optimize();
    }

}

GraphWrapperG2O::~GraphWrapperG2O()
{
    delete _so;
}


g2o::SparseOptimizer *GraphWrapperG2O::newSparseOptimizer(bool verbose) const
{
    g2o::SparseOptimizer *so = new g2o::SparseOptimizer;
    SlamLinearSolver *linearSolver = new SlamLinearSolver;
    linearSolver->setBlockOrdering(false);
    SlamBlockSolver *blockSolver = new SlamBlockSolver(linearSolver);

    so->setAlgorithm(new SlamOptimizationAlgorithm(blockSolver));
    so->setVerbose(verbose);
    return so;
}

GraphWrapperISAM *GraphWrapperG2O::toISAM()
{
    GraphWrapperISAM *isam = new GraphWrapperISAM;

    for(auto v: _so->indexMapping()) {
        isam->addVertex(v->id(), _vertexLookup[v]->estimate());
    }

    for(auto e: _so->edges()) {
        Eigen::MatrixXd info;
#ifndef G2S_QUATERNIONS
        info = _edgeLookup[e]->information();
#else /* G2S_QUATERNIONS */
        if(_edgeLookup[e]->measurement().is2d()) {
            info = _edgeLookup[e]->information();
        } else {
            info = EdgeSE3ISAM::eulerInformationFromQuat(_edgeLookup[e]->information());
        }
#endif /* G2S_QUATERNIONS */
        if(e->vertex(0)->id() == 0) {
            isam->addEdge(e->vertex(1)->id(),
                    _edgeLookup[e]->measurement(), info);
        } else {
            isam->addEdge(e->vertex(0)->id(), e->vertex(1)->id(),
                    _edgeLookup[e]->measurement(), info);
        }
    }

    isam->optimize();

    return isam;
}

GraphWrapper::Vertex *GraphWrapperG2O::vertex(int id)
{
    std::vector<GraphWrapper::Vertex *>::iterator it =
            std::lower_bound(_vertices.begin(), _vertices.end(), id, vertexCompare);
    return (it != _vertices.end() && (*it)->id() == id) ? *it : NULL;
}

void GraphWrapperG2O::addVertex(int id, const IsometryXd &init)
{
    g2o::OptimizableGraph::Vertex *nv;
    if(init.is2d()) {
        nv = new VertexSE2Type;
        static_cast<VertexSE2Type *>(nv)->setEstimate(init.se2());
    } else {
        nv = new VertexSE3Type;
        static_cast<VertexSE3Type *>(nv)->setEstimate(init.se3());
    }
    nv->setId(id);
    _so->addVertex(nv);
    _vertices.push_back(new Vertex(nv, _edgeLookup));
    _vertexLookup[nv] = _vertices.back();
}

void GraphWrapperG2O::addEdge(int from, int to, const IsometryXd &meas, const Eigen::MatrixXd &info)
{
    g2o::OptimizableGraph::Edge *ne;
    if(meas.is2d()) {
        ne = new EdgeSE2Type;
        static_cast<EdgeSE2Type *>(ne)->setMeasurement(meas.se2());
        static_cast<EdgeSE2Type *>(ne)->setInformation(info);
    }
    else {
        ne = new EdgeSE3Type;
        static_cast<EdgeSE3Type *>(ne)->setMeasurement(meas.se3());
        static_cast<EdgeSE3Type *>(ne)->setInformation(info);
    } 
    ne->setVertex(0, _so->vertex(from));
    ne->setVertex(1, _so->vertex(to));
    _so->addEdge(ne);
    _edgeLookup[ne] = new Edge(ne, _vertexLookup);
}


void GraphWrapperG2O::optimize()
{
    _so->vertex(0)->setFixed(true);
    if(_robust) {
        for(auto edge: _so->edges()) {
            g2o::OptimizableGraph::Edge *e = static_cast<g2o::OptimizableGraph::Edge *>(edge);
            g2o::RobustKernelCauchy *cauchy = new g2o::RobustKernelCauchy;
            cauchy->setDelta(1);
            e->setRobustKernel(cauchy);
        }
        _so->initializeOptimization();
        _so->optimize(50);
        for(auto edge: _so->edges()) {
            g2o::OptimizableGraph::Edge *e = static_cast<g2o::OptimizableGraph::Edge *>(edge);
            e->setRobustKernel(NULL);
        }
    }
    _so->initializeOptimization();
    _so->optimize(50);
}

void GraphWrapperG2O::optimizeFromId(int first_vertex_id)
{
    _so->vertex(first_vertex_id)->setFixed(true);
    if(_robust) {
        for(auto edge: _so->edges()) {
            g2o::OptimizableGraph::Edge *e = static_cast<g2o::OptimizableGraph::Edge *>(edge);
            g2o::RobustKernelCauchy *cauchy = new g2o::RobustKernelCauchy;
            cauchy->setDelta(1);
            e->setRobustKernel(cauchy);
        }
        _so->initializeOptimization();
        _so->optimize(50);
        for(auto edge: _so->edges()) {
            g2o::OptimizableGraph::Edge *e = static_cast<g2o::OptimizableGraph::Edge *>(edge);
            e->setRobustKernel(NULL);
        }
    }
    _so->initializeOptimization();
    _so->optimize(50);
}

void GraphWrapperG2O::push()
{
    _so->push();
}

void GraphWrapperG2O::pop()
{
    _so->pop();
}

void GraphWrapperG2O::optimizeWithReference(GraphWrapperG2O *other)
{
    _so->vertex(0)->setFixed(true);
    //_so->setVerbose(true);
    for(auto v: other->vertices()) {
        if(_so->vertex(v->id())) {
            _so->vertex(v->id())->setFixed(true);
            setEstimate(v->id(), v->estimate());
        }
    }
    _so->initializeOptimization();
    if(_so->indexMapping().size() > 0)
        _so->optimize(50);
    
    for(auto v: _vertices) {
        if(_so->vertex(v->id())) {
            _so->vertex(v->id())->setFixed(false);
        }
    }   
    _so->vertex(0)->setFixed(true);
    
    _so->initializeOptimization();
    _so->optimize(50);
    //_so->setVerbose(false);
}

static bool vertexIdCompare(
        const GraphWrapper::Vertex *v1,
        const GraphWrapper::Vertex *v2) {
    return v1->id() < v2->id();
}

GraphWrapperG2O *GraphWrapperG2O::clonePortion(int maxid)
{
    GraphWrapperG2O *gw = new GraphWrapperG2O(_verbose, _useGLC);

    for(auto v: _so->vertices()) {
        if(v.second->id() <= maxid) {
            gw->addVertex(v.second->id(), _vertexLookup[v.second]->estimate());
        }
    }

    for(auto e: _so->edges()) {
        if(e->vertex(0)->id() <= maxid && e->vertex(1)->id() <= maxid) {
            gw->addEdge(
                    e->vertex(0)->id(), e->vertex(1)->id(),
                    _edgeLookup[e]->measurement(), _edgeLookup[e]->information());
        }
    }

    std::sort(gw->_vertices.begin(), gw->_vertices.end(), vertexIdCompare);

    gw->optimize();
    return gw;
}

Eigen::VectorXd GraphWrapperG2O::estimate()
{
    return stack(_so);
}

Eigen::MatrixXd GraphWrapperG2O::information()
{
    return sparseInformation();
}

Eigen::MatrixXd GraphWrapperG2O::covariance()
{
    int nvars = fullSize(_so);
    Eigen::MatrixXd info = information();
    return info.llt().solve(Eigen::MatrixXd::Identity(nvars, nvars));
}

static bool vertexLessThan(
        const g2o::OptimizableGraph::Vertex *v1,
        const g2o::OptimizableGraph::Vertex *v2)
{
    return v1->id() < v2->id();
}

Eigen::SparseMatrix<double> GraphWrapperG2O::sparseInformation() const
{
    Eigen::SparseMatrix<double> info = sparse(static_cast<SlamBlockSolver *>(
            static_cast<g2o::OptimizationAlgorithmWithHessian *>(
                    _so->solver())->solver())->information());
    g2o::OptimizableGraph::VertexContainer verts = _so->indexMapping();

    // TODO: fix this
    for(size_t i = 0; i < verts.size() - 1; i++) {
        if(verts[i]->id() > verts[i + 1]->id())
            std::cerr << "NOPE! No good" << std::endl;
    }

    return info;
}

void GraphWrapperG2O::marginalizeNoOptimize(const std::vector<int> &which, const SparsityOptions &options)
{
    std::vector<g2o::OptimizableGraph::Vertex *> removeList;

    for(int id: which) {
        /* Search for vertex object in array */
        std::vector<GraphWrapper::Vertex *>::iterator it =
                std::lower_bound(_vertices.begin(), _vertices.end(), id, vertexCompare);
        assert(it != _vertices.end() && (*it)->id() == id &&
                "vertex needs to exist in order to be marginalized");

        /* Remove vertex object from map and array (but not from g2o!) */
        g2o::OptimizableGraph::Vertex *toRemove = _so->vertex(id);
        removeList.push_back(toRemove);
        _vertices.erase(it);
        _vertexLookup.erase(toRemove);

        /**
         * Canyu Le
         * To fix memory leak, edge objects will be removed in "vertex_remover.cpp" VertexRemover::updateInputGraph() function.
         * */
        // -----------------------
        /* Remove all edge objects attached to the vertex to marginalize (but not from g2o!) ERROR: below code will cause inconsistent edge remove */
//        std::cout<<"remove dummy edge\n";
//        for(g2o::HyperGraph::Edge *e: toRemove->edges()) {
//            for(auto v:e->vertices())
//                std::cout<<v->id()<<", ";
//            std::cout<<"\n";
//            _edgeLookup.erase(e);
//        }
        // -----------------------
    }

    VertexRemover vr;
    if(_useGLC) {
        vr.registerTopologyProvider(new TopologyProviderGLC);
    } else {
        vr.registerTopologyProvider(new TopologyProviderSE2);
        vr.registerTopologyProvider(new TopologyProviderSE2ISAM);
        vr.registerTopologyProvider(new TopologyProviderSE3);
        vr.registerTopologyProvider(new TopologyProviderSE3ISAM);
    }
    vr.setGraph(_so);
    vr.setEdgeMap(&_edgeLookup);
    vr.setSparsityOptions(options);

//    for(auto v: removeList) {

    g2o::OptimizableGraph::EdgeContainer added = vr.remove(removeList);

    for(auto e: added) {
        _edgeLookup[e] = new Edge(e, _vertexLookup);
    }

//    }
}

void GraphWrapperG2O::marginalize(const std::vector<int> &which, const SparsityOptions &options)
{
    marginalizeNoOptimize(which, options);
//    for(int v: which) {
//        std::vector<int> toRemove(1, v);
//        marginalizeNoOptimize(toRemove, options);
//    }
    optimize();
}

#include <g2o/core/factory.h>

void GraphWrapperG2O::write(std::ofstream &s)
{
    _so->save(s);
}

void GraphWrapperG2O::computeIndices(
        GraphWrapper *other,
        std::vector<int> &indicesKeep,
        std::vector<int> &indicesMarginalize)
{
    const std::vector<GraphWrapper::Vertex *> &vertices = other->vertices();
    indicesKeep.clear();
    indicesMarginalize.clear();

    for(int i = 1, j = 0, k = 0; i < _vertices.size(); i++) { /* Auto-skip first node */
        int id = _vertices[i]->id();
        if(j < vertices.size() && vertices[j]->id() == 0) j++;

        /* int on X::Dimension to avoid int& reference. Without it g++ debug mode doesn't link */
        int dim = _vertices[i]->is2d() ?
                int(EdgeSE2Type::Dimension) : int(EdgeSE3Type::Dimension);

        if(j == vertices.size() || id < vertices[j]->id()) {
            for(int t = 0; t < dim; t++)
                indicesMarginalize.push_back(k + t);
        } else {
            for(int t = 0; t < dim; t++)
                indicesKeep.push_back(k + t);
            j++;
        }
        k += dim;
    }
}

std::ostream &operator<<(std::ostream &s, GraphWrapper::Vertex *v)
{
    s << v->id();
}

double GraphWrapperG2O::chi2(GraphWrapper *other)
{
    _so->push();

    for(GraphWrapper::Vertex *v: other->vertices()) {
        setEstimate(v->id(), v->estimate());
        _so->vertex(v->id())->setFixed(true);
    }

    optimize();

    double chi2val = _so->chi2();

    for(GraphWrapper::Vertex *v: other->vertices()) {
        _so->vertex(v->id())->setFixed(false);
    }

    _so->vertex(0)->setFixed(true);

    _so->pop();

    return chi2val;
}


double GraphWrapperG2O::kullbackLeibler(GraphWrapper *other)
{
    Eigen::VectorXd fullMean = estimate();
    Eigen::SparseMatrix<double> fullInfo = sparseInformation();
    std::vector<int> indicesKeep, indicesMarginalize;

    computeIndices(other, indicesKeep, indicesMarginalize);

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > chol(selectVariables(fullInfo, indicesMarginalize));
    Eigen::SparseMatrix<double> mixed = selectVariables(fullInfo, indicesMarginalize, indicesKeep);
    Eigen::MatrixXd selectedInfo = selectVariables(fullInfo, indicesKeep);
    selectedInfo -= mixed.transpose() * chol.solve(mixed);

    double kld = kullbackLeiblerDivergence(
            estimateDifference(other), other->information(), selectedInfo, InformationInformation);

    return kld;
}

Eigen::VectorXd GraphWrapperG2O::estimateDifference(GraphWrapper *other)
{
    const std::vector<GraphWrapper::Vertex *> &vertices = other->vertices();
    Eigen::VectorXd diff(fullSize(static_cast<GraphWrapperG2O *>(other)->_so));
    for(int i = 1, j = 0, k = 0; i < _vertices.size(); i++) { /* Auto-skip first node */
        if(j < vertices.size() && vertices[j]->id() == 0) j++;
        if(j == vertices.size()) {
            break;
        } else if(_vertices[i]->id() == vertices[j]->id()) {
            if(_vertices[i]->is2d()) {
                Eigen::Vector3d d = _vertices[i]->estimate().vector() -
                        vertices[j]->estimate().vector();
                d[2] = g2o::normalize_theta(d[2]);
                diff.segment<3>(k) = d;
                k += 3;
            } else {
                IsometryXd isod = _vertices[i]->estimate().inverse() * vertices[j]->estimate();
                Eigen::VectorXd d = isod.vector(IsometryXd::CondensedQuaternion);
                diff.segment<6>(k) = d;
                k += 6;
            }
            j++;
        }
    }
    return diff;
}

void GraphWrapperG2O::debugPrint(std::ostream &s) const
{
    s << "+ vertices: ";
    for(auto v: _so->vertices()) {
        s << v.first << " ";
    }
    s << "\n+ edges: ";
    for(auto e: _so->edges()) {
        auto vertices = e->vertices();
        s << "(";
        for(size_t i = 0; i < vertices.size(); i++) {
            s << vertices[i]->id();
            if(i < vertices.size() - 1) s << ",";
        }
        s << ") ";
    }
    s << std::endl;
}

void GraphWrapperG2O::saveInformation(const char *fname)
{
    std::ofstream f(fname);
    std::vector<Eigen::Triplet<double> > trs = triplets(sparseInformation());
    for(auto &t: trs) {
        f << t.row() + 1 << "\t" << t.col() + 1 << "\t" << t.value() << std::endl;
    }
    f.close();
}

void GraphWrapperG2O::printStats(std::ostream &s) const
{
    Eigen::SparseMatrix<double> info = sparseInformation();
    double fillin = info.nonZeros() / double(info.rows() * info.cols());
    s << "nodes = " << _so->vertices().size() - 1 << "; edges = " << _so->edges().size() <<
            "; fillin = " << fillin * 100 << "%";
}

void GraphWrapperG2O::setEstimate(int vertexid, const IsometryXd &est)
{
    if(est.is2d()) {
        static_cast<VertexSE2Type *>(_so->vertex(vertexid))->setEstimate(est.se2());
    } else {
        static_cast<VertexSE3Type *>(_so->vertex(vertexid))->setEstimate(est.se3());
    }
}
