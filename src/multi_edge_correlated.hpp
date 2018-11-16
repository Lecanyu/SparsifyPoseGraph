/*
 * multi_edge_correlated.hpp
 *
 *  Created on: Oct 15, 2014
 *      Author: mazuran
 */

#ifndef MULTI_EDGE_CORRELATED_HPP_
#define MULTI_EDGE_CORRELATED_HPP_

#include <sstream>

template <typename EdgeType>
MultiEdgeCorrelated<EdgeType>::MultiEdgeCorrelated() :
    Base(), _mockedge(new EdgeType)
{
    resize(0);
    _mockedge->setInformation(Eigen::MatrixXd::Identity(
            EdgeType::Dimension, EdgeType::Dimension));
}

template <typename EdgeType>
MultiEdgeCorrelated<EdgeType>::~MultiEdgeCorrelated()
{
    delete _mockedge;
}

template <typename EdgeType>
void MultiEdgeCorrelated<EdgeType>::addMeasurement(
            const g2o::HyperGraph::VertexContainer &verts,
            const Measurement &meas)
{
    std::list<g2o::HyperGraph::Vertex *> toAdd;
    std::list<int> indices;
    std::vector<int> toAddIndices;
    toAddIndices.reserve(verts.size());
    int counter = _vertices.size();
    for(g2o::HyperGraph::VertexContainer::const_iterator it = verts.begin();
            it != verts.end(); ++it) {
        g2o::HyperGraph::VertexContainer::const_iterator where = std::find(
                _vertices.begin(), _vertices.end(), *it);
        if(where == _vertices.end()) {
            toAdd.push_back(*it);
            toAddIndices.push_back(counter);
            indices.push_back(counter++);
        } else {
            indices.push_back(where - _vertices.begin());
        }
    }

    if(toAdd.size() > 0) {
        resize(_vertices.size() + toAdd.size());
        int i = 0;
        for(std::list<g2o::HyperGraph::Vertex *>::const_iterator it = toAdd.begin();
                it != toAdd.end(); ++it, i++) {
            _vertices[toAddIndices[i]] = *it;
        }
    }

    _mappings.push_back(indices);
    _measurement.push_back(meas);
}

template <typename EdgeType>
void MultiEdgeCorrelated<EdgeType>::computeError()
{
    int k = 0;
    MappingsType::const_iterator itv = _mappings.begin();
    for(typename std::list<Measurement>::const_iterator itm = _measurement.begin();
            itm != _measurement.end() && itv != _mappings.end(); ++itm, ++itv) {
        _mockedge->vertices() = toVertexContainer(*itv);
        _mockedge->setMeasurement(*itm);
        _mockedge->computeError();
        _error.template segment<EdgeType::Dimension>(k) = _mockedge->error();
        k += EdgeType::Dimension;
    }
}

template <typename EdgeType>
bool MultiEdgeCorrelated<EdgeType>::setMeasurementFromState()
{
    int k = 0;
    MappingsType::const_iterator itv = _mappings.begin();
    for(typename std::list<Measurement>::iterator itm = _measurement.begin();
            itm != _measurement.end() && itv != _mappings.end(); ++itm, ++itv) {
        _mockedge->vertices() = toVertexContainer(*itv);
        if(!_mockedge->setMeasurementFromState()) {
            return false;
        }
        *itm = _mockedge->measurement();
    }

    return true;
}

template <typename EdgeType>
void MultiEdgeCorrelated<EdgeType>::linearizeOplus()
{
    size_t nedgeverts = 0;
    g2o::JacobianWorkspace jw;
    int maxVertexSize = 0, numVertices = 0;
   
    for(MappingsType::const_iterator itv = _mappings.begin();
            itv != _mappings.end(); ++itv) {
        numVertices = std::max(numVertices, (int) itv->size());
    }

    for(g2o::HyperGraph::Vertex *v: _vertices) {
        maxVertexSize = std::max(static_cast<g2o::OptimizableGraph::Vertex *>(
                v)->dimension(), maxVertexSize);
    }
    
    jw.updateSize(numVertices, EdgeType::Dimension * maxVertexSize);
    jw.allocate();

    int i = 0, k = 0;
    for(g2o::HyperGraph::VertexContainer::const_iterator it = _vertices.begin();
            it != _vertices.end(); ++it, i++) {
        _jacobianOplus[i] = Eigen::MatrixXd::Zero(
                _measurement.size() * EdgeType::Dimension,
                static_cast<g2o::OptimizableGraph::Vertex *>(*it)->dimension());
    }

    MappingsType::const_iterator itv = _mappings.begin();
    for(typename std::list<Measurement>::const_iterator itm = _measurement.begin();
            itm != _measurement.end() && itv != _mappings.end();
            ++itm, ++itv, k += EdgeType::Dimension) {
        _mockedge->vertices() = toVertexContainer(*itv);
        _mockedge->setMeasurement(*itm);
        static_cast<g2o::OptimizableGraph::Edge *>(_mockedge)->linearizeOplus(jw);

        i = 0;
        for(std::list<int>::const_iterator it = itv->begin();
                it != itv->end(); ++it, i++) {
            int vdim = static_cast<g2o::OptimizableGraph::Vertex *>(_vertices[*it])->dimension();
            Eigen::Map<Eigen::MatrixXd> Jv(jw.workspaceForVertex(i), EdgeType::Dimension, vdim);
            _jacobianOplus[*it].block(k, 0, EdgeType::Dimension, vdim) = Jv;
        }
    }
}

template <typename EdgeType>
void MultiEdgeCorrelated<EdgeType>::initialEstimate(
        const g2o::HyperGraph::VertexSet &fixed,
        g2o::OptimizableGraph::Vertex *toEstimate)
{
    MappingsType::const_iterator itv = _mappings.begin();
    for(typename std::list<Measurement>::const_iterator itm = _measurement.begin();
            itm != _measurement.end() && itv != _mappings.end(); ++itm, ++itv) {
        g2o::HyperGraph::VertexContainer vc = toVertexContainer(*itv);
        if(setMembership(vc, toEstimate)) {
            _mockedge->vertices() = vc;
            _mockedge->setMeasurement(*itm);
            _mockedge->initialEstimate(setIntersection(fixed, vc), toEstimate);
        }
    }
}

template <typename EdgeType>
double MultiEdgeCorrelated<EdgeType>::initialEstimatePossible(
        const g2o::HyperGraph::VertexSet &fixed,
        g2o::OptimizableGraph::Vertex *toEstimate)
{
    double cost = 0;
    MappingsType::const_iterator itv = _mappings.begin();
    for(typename std::list<Measurement>::const_iterator itm = _measurement.begin();
            itm != _measurement.end() && itv != _mappings.end(); ++itm, ++itv) {
        g2o::HyperGraph::VertexContainer vc = toVertexContainer(*itv);
        if(setMembership(vc, toEstimate)) {
            _mockedge->vertices() = vc;
            _mockedge->setMeasurement(*itm);
            double thiscost = _mockedge->initialEstimatePossible(
                    setIntersection(fixed, vc), toEstimate);
            if(thiscost < 0)
                return -1;
            cost += thiscost;
        }
    }
    return cost;
}

template <typename EdgeType>
bool MultiEdgeCorrelated<EdgeType>::read(std::istream& is)
{
    int nmeas = 0, nrelevant = 0, n = _mockedge->dimension();
    is >> nmeas >> nrelevant;
    setMeasurementCount(nmeas);

    for(int i = 0; i < nmeas; i++) {
        /*
         * TODO: This is a hack. It adds the information matrix to the string
         * to be read, however g2o doesn't enforce read() to load n*(n+1)/2
         * entries for the information matrix at the end, hence it may not work.
         * But it should work for all default g2o types.
         */
        std::stringstream ss;
        for(int j = 0; j < nrelevant; j++) {
            std::string str;
            is >> str;
            ss << str << " ";
        }
        for(int j = 0; j < n; j++) {
            ss << "1 ";
            for(int k = j + 1; k < n; k++) {
                ss << "0 ";
            }
        }
        ss.seekg(0);

        if(!_mockedge->read(ss)) return false;
        _measurement.push_back(_mockedge->measurement());
    }

    for(int i = 0; i < _information.rows(); i++) {
        for(int j = i; j < _information.cols(); j++) {
            is >> _information(i, j);
        }
    }

    _information.template triangularView<Eigen::StrictlyLower>() =
            _information.template triangularView<Eigen::StrictlyUpper>().transpose();

    return true;
}

template <typename EdgeType>
bool MultiEdgeCorrelated<EdgeType>::write(std::ostream& os) const
{
    /* If _mockedge->write() returns false, this avoids writing to os */
    std::stringstream everything;
    everything << "|| " << _measurement.size();
    int n = _mockedge->dimension();

    for(typename std::list<Measurement>::const_iterator it = _measurement.begin();
            it != _measurement.end(); ++it) {
        _mockedge->setMeasurement(*it);
        std::stringstream ss;
        if(!_mockedge->write(ss)) return false;

        /*
         * TODO: This is a hack. It removes the information matrix from the
         * written string, however g2o doesn't enforce write() to append n*(n+1)/2
         * entries for the information matrix at the end, hence it may not work.
         * But it should work for all default g2o types.
         */
        ss.seekg(0);
        std::vector<std::string> tokens;
        while(ss.good()) {
            std::string str;
            ss >> str;
            tokens.push_back(str);
        }
        int nrelevant = tokens.size() - n * (n + 1) / 2;
        if(it == _measurement.begin())     everything << " " << nrelevant;
        for(int i = 0; i < nrelevant; i++) everything << " " << tokens[i];
    }

    os << everything.str();

    for(int i = 0; i < _information.rows(); i++) {
        for(int j = i; j < _information.cols(); j++) {
            os << " " << _information(i, j);
        }
    }

    return true;
}

template <typename EdgeType>
g2o::HyperGraph::VertexSet MultiEdgeCorrelated<EdgeType>::setIntersection(
            const g2o::HyperGraph::VertexSet &a,
            const g2o::HyperGraph::VertexContainer &b)
{
    g2o::HyperGraph::VertexSet ret;
    for(g2o::HyperGraph::VertexContainer::const_iterator it = b.begin();
            it != b.end(); ++it) {
        if(a.count(*it) > 0) ret.insert(*it);
    }
    return ret;
}

template <typename EdgeType>
bool MultiEdgeCorrelated<EdgeType>::setMembership(
            const g2o::HyperGraph::VertexContainer &a,
            const g2o::HyperGraph::Vertex *v)
{
    return std::find(a.begin(), a.end(), v) != a.end();
}

template <typename EdgeType>
g2o::HyperGraph::VertexContainer MultiEdgeCorrelated<EdgeType>::toVertexContainer(
            const std::list<int> &vmap) const
{
    g2o::HyperGraph::VertexContainer vs;
    for(std::list<int>::const_iterator it = vmap.begin();
            it != vmap.end(); ++it) {
        vs.push_back(_vertices[*it]);
    }
    return vs;
}


#endif /* MULTI_EDGE_CORRELATED_HPP_ */
