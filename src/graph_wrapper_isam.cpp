/*
 * graph_wrapper_isam.cpp
 *
 *  Created on: 15/dic/2013
 *      Author: Mladen Mazuran
 */

#include "graph_wrapper_isam.h"
#include "utils.h"
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <isam/slam2d.h>
#include <isam/slam3d.h>
#include <isam/glc.h>

/* ---------------------------------------------------------------------------------------------- */
/*                                    GraphWrapperISAM::Vertex                                    */
/* ---------------------------------------------------------------------------------------------- */

GraphWrapperISAM::Vertex::Vertex(int id, isam::Node *v, const EdgeMap &l) :
        _id(id), _v(v), _l(&l)
{
}

std::vector<const GraphWrapper::Edge *> GraphWrapperISAM::Vertex::edges() const
{
    std::vector<const GraphWrapper::Edge *> ret;
    for(const isam::Factor *e: _v->factors()) {
        ret.push_back(_l->at(e));
    }
    return ret;
}

IsometryXd GraphWrapperISAM::Vertex::estimate() const
{
    const isam::Pose2d_Node *vse2 = dynamic_cast<const isam::Pose2d_Node *>(_v);
    const isam::Pose3d_Node *vse3 = dynamic_cast<const isam::Pose3d_Node *>(_v);
    assert((vse2 || vse3) && "Can handle only SE2 or SE3 vertices for now");
    return vse2 ? IsometryXd(vse2->value()) : IsometryXd(vse3->value());
}

bool GraphWrapperISAM::Vertex::is2d() const
{
    return dynamic_cast<isam::Pose2d_Node *>(_v);
}

/* ---------------------------------------------------------------------------------------------- */
/*                                     GraphWrapperISAM::Edge                                     */
/* ---------------------------------------------------------------------------------------------- */

GraphWrapperISAM::Edge::Edge(isam::Factor *e, const GraphWrapperISAM::VertexMap &l) :
        _e(e), _l(&l)
{
}

std::vector<const GraphWrapper::Vertex *> GraphWrapperISAM::Edge::vertices() const
{
    std::vector<const GraphWrapper::Vertex *> ret;
    for(const isam::Node *v: _e->nodes()) {
        ret.push_back(_l->at(v));
    }
    return ret;
}

IsometryXd GraphWrapperISAM::Edge::measurement() const
{
    const isam::Pose2d_Pose2d_Factor *ese2 = dynamic_cast<const isam::Pose2d_Pose2d_Factor *>(_e);
    const isam::Pose3d_Pose3d_Factor *ese3 = dynamic_cast<const isam::Pose3d_Pose3d_Factor *>(_e);
    assert((ese2 || ese3) && "Can handle only SE2 or SE3 binary edges for now");
    return ese2 ? IsometryXd(ese2->measurement()) : IsometryXd(ese3->measurement());
}

Eigen::MatrixXd GraphWrapperISAM::Edge::information() const
{
    const isam::Pose2d_Pose2d_Factor *ese2 = dynamic_cast<const isam::Pose2d_Pose2d_Factor *>(_e);
    const isam::Pose3d_Pose3d_Factor *ese3 = dynamic_cast<const isam::Pose3d_Pose3d_Factor *>(_e);
    assert((ese2 || ese3) && "Can handle only SE2 or SE3 binary edges for now");
    Eigen::MatrixXd sqrtinfo = ese2 ?
            Eigen::MatrixXd(ese2->sqrtinf()) : Eigen::MatrixXd(ese3->sqrtinf());
    Eigen::MatrixXd info = sqrtinfo.transpose() * sqrtinfo;
    return 0.5 * (info + info.transpose());
}

/* ---------------------------------------------------------------------------------------------- */
/*                                        GraphWrapperISAM                                        */
/* ---------------------------------------------------------------------------------------------- */

GraphWrapperISAM::GraphWrapperISAM() :
    _slam(new isam::Slam)
{
    isam::Properties prop;
#ifdef G2S_VERBOSE_DEBUG
    prop.verbose = true;
    prop.quiet = false;
#else /* G2S_VERBOSE_DEBUG */
    prop.verbose = false;
    prop.quiet = true;
#endif /* G2S_VERBOSE_DEBUG */

    prop.mod_solve = isam::LEVENBERG_MARQUARDT;
    prop.epsilon_abs = 1e-8; //1e-5;
    prop.epsilon_rel = 1e-8; //1e-6;
    prop.mod_batch = 1;
    _slam->set_properties(prop);
}

GraphWrapperISAM::~GraphWrapperISAM()
{
    /* Delete isam optimizer, nodes, and factors */
    for(auto f: _slam->get_factors()) {
        delete f;
    }
    for(auto n: _slam->get_nodes()) {
        delete n;
    }
    delete _slam;
}

void GraphWrapperISAM::addVertex(int id, const IsometryXd &init)
{
    isam::Node *n;
    if(init.is2d()) {
        n = new isam::Pose2d_Node;
        static_cast<isam::Pose2d_Node *>(n)->init(isam::Pose2d(init.vector()));
    } else {
        n = new isam::Pose3d_Node;
        static_cast<isam::Pose3d_Node *>(n)->init(isam::Pose3d(init.se3()));
    }
    _slam->add_node(n);
    _forwardLookup[id] = n;
    _reverseLookup[n]  = id;
    _vertices.push_back(new Vertex(id, n, _edgeLookup));
    _vertexLookup[n] = _vertices.back();
}

void GraphWrapperISAM::addEdge(int from, int to, const IsometryXd &meas, const Eigen::MatrixXd &info)
{
    isam::Factor *f;
    if(meas.is2d()) {
        f = new isam::Pose2d_Pose2d_Factor(
            static_cast<isam::Pose2d_Node *>(_forwardLookup[from]),
            static_cast<isam::Pose2d_Node *>(_forwardLookup[to]),
            isam::Pose2d(meas.vector()), isam::Information(info));
    } else {
        f = new isam::Pose3d_Pose3d_Factor(
            static_cast<isam::Pose3d_Node *>(_forwardLookup[from]),
            static_cast<isam::Pose3d_Node *>(_forwardLookup[to]),
            isam::Pose3d(meas.se3()), isam::Information(info));
    }
    _slam->add_factor(f);
    _edgeLookup[f] = new Edge(f, _vertexLookup);
}

void GraphWrapperISAM::addEdge(int to, const IsometryXd &meas, const Eigen::MatrixXd &info)
{
    isam::Factor *f;
    if(meas.is2d()) {
        f = new isam::Pose2d_Factor(
            static_cast<isam::Pose2d_Node *>(_forwardLookup[to]),
            isam::Pose2d(meas.vector()), isam::Information(info));
    } else {
        f = new isam::Pose3d_Factor(
            static_cast<isam::Pose3d_Node *>(_forwardLookup[to]),
            isam::Pose3d(meas.se3()), isam::Information(info));
    }
    _slam->add_factor(f);
    _edgeLookup[f] = new Edge(f, _vertexLookup);
}

void GraphWrapperISAM::addEdgeSqrt(int from, int to, const IsometryXd &meas, const Eigen::MatrixXd &sqrtinfo)
{
    isam::Factor *f;
    if(meas.is2d()) {
        f = new isam::Pose2d_Pose2d_Factor(
            static_cast<isam::Pose2d_Node *>(_forwardLookup[from]),
            static_cast<isam::Pose2d_Node *>(_forwardLookup[to]),
            isam::Pose2d(meas.vector()), isam::SqrtInformation(sqrtinfo));
    } else {
        f = new isam::Pose3d_Pose3d_Factor(
            static_cast<isam::Pose3d_Node *>(_forwardLookup[from]),
            static_cast<isam::Pose3d_Node *>(_forwardLookup[to]),
            isam::Pose3d(meas.se3()), isam::SqrtInformation(sqrtinfo));
    }
    _slam->add_factor(f);
    _edgeLookup[f] = new Edge(f, _vertexLookup);
}

void GraphWrapperISAM::addEdgeSqrt(int to, const IsometryXd &meas, const Eigen::MatrixXd &sqrtinfo)
{
    isam::Factor *f;
    if(meas.is2d()) {
        f = new isam::Pose2d_Factor(
            static_cast<isam::Pose2d_Node *>(_forwardLookup[to]),
            isam::Pose2d(meas.vector()), isam::SqrtInformation(sqrtinfo));
    } else {
        f = new isam::Pose3d_Factor(
            static_cast<isam::Pose3d_Node *>(_forwardLookup[to]),
            isam::Pose3d(meas.se3()), isam::SqrtInformation(sqrtinfo));
    }
    _slam->add_factor(f);
    _edgeLookup[f] = new Edge(f, _vertexLookup);
}

void GraphWrapperISAM::optimize()
{
    _slam->batch_optimization();
#ifdef G2S_VERBOSE_DEBUG
    _slam->print_stats();
#endif /* G2S_VERBOSE_DEBUG */
}

GraphWrapperISAM *GraphWrapperISAM::clonePortion(int maxid)
{
    GraphWrapperISAM *gw = new GraphWrapperISAM;
    for(auto pair: _forwardLookup) {
        if(pair.first <= maxid) {
            gw->addVertex(pair.first, _vertexLookup[pair.second]->estimate());
        }
    }

    for(auto e: _slam->get_factors()) {
        if(e->nodes().size() == 1) {
            isam::Pose2d_Factor *ese2 = dynamic_cast<isam::Pose2d_Factor *>(e);
            isam::Pose3d_Factor *ese3 = dynamic_cast<isam::Pose3d_Factor *>(e);
            if(ese2) {
                gw->addEdgeSqrt(_reverseLookup[e->nodes()[0]], IsometryXd(ese2->measurement()), ese2->sqrtinf());
            } else {
                gw->addEdgeSqrt(_reverseLookup[e->nodes()[0]], IsometryXd(ese3->measurement()), ese3->sqrtinf());
            }
        } else if(_reverseLookup[e->nodes()[0]] <= maxid && _reverseLookup[e->nodes()[1]] <= maxid) {
            isam::Pose2d_Pose2d_Factor *ese2 = dynamic_cast<isam::Pose2d_Pose2d_Factor *>(e);
            isam::Pose3d_Pose3d_Factor *ese3 = dynamic_cast<isam::Pose3d_Pose3d_Factor *>(e);
            if(ese2) {
                gw->addEdgeSqrt(_reverseLookup[e->nodes()[0]], _reverseLookup[e->nodes()[1]], IsometryXd(ese2->measurement()), ese2->sqrtinf());
            } else {
                gw->addEdgeSqrt(_reverseLookup[e->nodes()[0]], _reverseLookup[e->nodes()[1]], IsometryXd(ese3->measurement()), ese3->sqrtinf());
            }
        }
    }

    gw->optimize();
    return gw;
}

Eigen::VectorXd GraphWrapperISAM::estimate()
{
    return stack(_slam);
}

Eigen::MatrixXd GraphWrapperISAM::information()
{
    //int nvars = fullSize(_slam);
    //return covariance().llt().solve(Eigen::MatrixXd::Identity(nvars, nvars));
    Eigen::SparseMatrix<double> R = sparse(_slam->get_R());
    Eigen::SparseMatrix<double> fullInfo = R.transpose() * R;
    return fullInfo;
}

Eigen::MatrixXd GraphWrapperISAM::covariance()
{
    return _slam->covariances().marginal(_slam->get_nodes());
}

static bool vertexCompare(const GraphWrapper::Vertex *v, int id)
{
    return v->id() < id;
}

void GraphWrapperISAM::marginalize(const std::vector<int> &which, const SparsityOptions &options)
{
    for(int id: which) {
        std::vector<GraphWrapper::Vertex *>::iterator it =
                std::lower_bound(_vertices.begin(), _vertices.end(), id, vertexCompare);
        assert(it != _vertices.end() && (*it)->id() == id &&
                "vertex needs to exist in order to be marginalized");

        /* Remove vertex object from map and array (but not from isam!) */
        isam::Node *toRemove = _forwardLookup[id];
        _vertices.erase(it);
        _vertexLookup.erase(toRemove);

        /* Remove all edge objects attached to the vertex to marginalize (but not from isam!) */
        for(isam::Factor *e: toRemove->factors()) {
            _edgeLookup.erase(e);
        }
        _reverseLookup.erase(toRemove);
        _forwardLookup.erase(id);

        std::vector<isam::Factor *> removed = glc_elim_factors(toRemove);

        /* TODO: the new is actually a memory leak, but it is not critical and I'm lazy, so for now
         * let's turn a blind eye */
        glc_remove_node(*_slam, toRemove, options.topology == SparsityOptions::Tree, new isam::GLC_RootShift);
    
        /* Now delete isam's objects */
        delete toRemove;
        for(auto factor: removed) {
            delete factor;
        }
    }
    optimize();
}

void GraphWrapperISAM::write(std::ofstream &s)
{
    _slam->write(s);
}

double GraphWrapperISAM::kullbackLeibler(GraphWrapper *other) {
    Eigen::VectorXd fullMean = estimate();
    Eigen::SparseMatrix<double> J = sparse(_slam->jacobian());
    Eigen::SparseMatrix<double> fullInfo = J.transpose() * J;

    auto vertices = other->vertices();
    std::vector<int> indicesKeep, indicesMarginalize;

    for(int i = 0, j = 0, k = 0; i < _vertices.size(); i++) {
        int id = _vertices[i]->id();
        if(j < vertices.size() && vertices[j]->id() == 0) j++;

        int dim = _vertices[i]->is2d() ? 3 : 6;
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

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > chol(selectVariables(fullInfo, indicesMarginalize));
    Eigen::SparseMatrix<double> mixed = selectVariables(fullInfo, indicesMarginalize, indicesKeep);
    Eigen::MatrixXd selectedInfo = selectVariables(fullInfo, indicesKeep);
    selectedInfo -= mixed.transpose() * chol.solve(mixed);

    double kld = kullbackLeiblerDivergence(
            estimateDifference(other), other->information(), selectedInfo, InformationInformation);

    return kld;
}

Eigen::VectorXd GraphWrapperISAM::estimateDifference(GraphWrapper *other)
{
    const std::vector<GraphWrapper::Vertex *> &vertices = other->vertices();
    Eigen::VectorXd diff(fullSize(static_cast<GraphWrapperISAM *>(other)->_slam));
    for(int i = 0, j = 0, k = 0; i < _vertices.size(); i++) {
        if(j < vertices.size() && vertices[j]->id() == 0) j++;
        if(j == vertices.size()) {
            break;
        } else if(_vertices[i]->id() == vertices[j]->id()) {
#ifdef G2S_MANIFOLD_ISAM
            int dim = vertices[j]->is2d() ? 3 : 6;
            IsometryXd isod = vertices[j]->estimate().inverse() * _vertices[i]->estimate();
            diff.segment(k, dim) = isod.vector(IsometryXd::EulerAnglesISAM);
            k += dim;
#else /* G2S_MANIFOLD_ISAM */
            if(_vertices[i]->is2d()) {
                Eigen::Vector3d d = _vertices[i]->estimate().vector() -
                        vertices[j]->estimate().vector();
                d[2] = g2o::normalize_theta(d[2]);
                diff.segment<3>(k) = d;
                k += 3;
            } else {
                /*
                Eigen::VectorXd d = _vertices[i]->estimate().vector(IsometryXd::EulerAnglesISAM) -
                        vertices[j]->estimate().vector(IsometryXd::EulerAnglesISAM);
                d[3] = g2o::normalize_theta(d[3]);
                d[4] = g2o::normalize_theta(d[4]);
                d[5] = g2o::normalize_theta(d[5]);
                */
                Eigen::Isometry3d iso = (_vertices[i]->estimate().inverse() * vertices[j]->estimate()).se3();
                Eigen::AngleAxisd aa(iso.rotation());
                Eigen::VectorXd d(6);
                d.segment<3>(0) = iso.translation();
                d.segment<3>(3) = aa.angle() * aa.axis();

                diff.segment<6>(k) = d;
                k += 6;
            }
#endif /* G2S_MANIFOLD_ISAM */
            j++;
        }
    }
    return diff;
}

void GraphWrapperISAM::saveInformation(const char *fname)
{
    std::ofstream f(fname);
    Eigen::SparseMatrix<double> R = sparse(_slam->get_R());
    Eigen::SparseMatrix<double> fullInfo = R.transpose() * R;
    std::vector<Eigen::Triplet<double> > trs = triplets(fullInfo);
    for(auto &t: trs) {
        f << t.row() + 1 << "\t" << t.col() + 1 << "\t" << t.value() << std::endl;
    }
    f.close();
}

void GraphWrapperISAM::debugPrint(std::ostream &s) const
{
    s << "+ vertices: ";
    for(auto n: _slam->get_nodes()) {
        s << _reverseLookup.at(n) << " ";
    }
    s << "\n+ edges: ";
    for(auto f: _slam->get_factors()) {
        std::vector<isam::Node *> nodes = f->nodes();
        s << "(";
        for(size_t i = 0; i < nodes.size(); i++) {
            s << _reverseLookup.at(nodes[i]);
            if(i < nodes.size() - 1) s << ",";
        }
        s << ") ";
    }
    s << std::endl;
}

void GraphWrapperISAM::printStats(std::ostream &s) const
{
    Eigen::SparseMatrix<double> J = sparse(_slam->jacobian());
    Eigen::SparseMatrix<double> info = J.transpose() * J;
    double fillin = info.nonZeros() / double(info.rows() * info.cols());
    s << "nodes = " << _slam->get_nodes().size() << "; edges = " << _slam->get_factors().size() <<
            "; fillin = " << fillin * 100 << "%";
}

void GraphWrapperISAM::setEstimate(int vertexid, const IsometryXd &est)
{
    if(est.is2d()) {
        static_cast<isam::Pose2d_Node *>(_forwardLookup[vertexid])->init(isam::Pose2d(est.vector()));
    } else {
        static_cast<isam::Pose3d_Node *>(_forwardLookup[vertexid])->init(isam::Pose3d(est.se3()));
    }
}
