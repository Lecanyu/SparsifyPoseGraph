/*
 * glc_edge.hpp
 *
 *  Created on: Oct 15, 2014
 *      Author: mazuran
 */

#include "glc_edge.h"
#include <sstream>
#include <g2o/core/factory.h>

GLCEdge::GLCEdge() :
    Base(), _reparam(NULL)
{
    resize(0);
}

GLCEdge::~GLCEdge()
{
    delete _reparam;
}

void GLCEdge::computeMeasurement()
{
    _measurement = _reparam->reparametrize(_vertices);
}

void GLCEdge::computeError()
{
    assert(_reparam && "needs GLCReparam");
    _error = _W * _reparam->reparametrize(_vertices, _measurement);
}

bool GLCEdge::setMeasurementFromState()
{
    return false;
}


void GLCEdge::linearizeOplus()
{
    assert(_reparam && "needs GLCReparam");
    g2o::MatrixXD J = _reparam->jacobian(_vertices, _measurement);
    for(size_t i = 0, k = 0, dim = 0; i < _vertices.size(); k += dim, i++) {
        dim = static_cast<g2o::OptimizableGraph::Vertex *>(
                _vertices[i])->dimension();
        _jacobianOplus[i] = _W * J.block(0, k, J.rows(), dim);
    }
}

void GLCEdge::initialEstimate(
        const g2o::HyperGraph::VertexSet &,
        g2o::OptimizableGraph::Vertex *)
{
}

double GLCEdge::initialEstimatePossible(
        const g2o::HyperGraph::VertexSet &,
        g2o::OptimizableGraph::Vertex *)
{
    return -1;
}

bool GLCEdge::read(std::istream& is)
{
    std::string reparamTag;
    is >> reparamTag;
    _reparam = static_cast<GLCReparam *>(g2o::Factory::instance()->construct(reparamTag));

    int errsize = 0, meassize = 0;
    is >> errsize >> meassize;
    setDimension(errsize, meassize);
    for(int i = 0; i < meassize; i++) {
        is >> _measurement[i];
    }

    for(int i = 0; i < errsize; i++) {
        for(int j = 0; j < meassize; j++) {
            is >> _W(i, j);
        }
    }

    for(int i = 0; i < _information.rows(); i++) {
        for(int j = i; j < _information.cols(); j++) {
            is >> _information(i, j);
        }
    }

    _information.triangularView<Eigen::StrictlyLower>() =
            _information.triangularView<Eigen::StrictlyUpper>().transpose();

    return true;
}

bool GLCEdge::write(std::ostream& os) const
{
    assert(_reparam && "needs GLCReparam");
    os << "|| " << g2o::Factory::instance()->tag(_reparam) << " ";

    os << _W.rows() << " " << _W.cols();

    for(int i = 0; i < _measurement.rows(); i++) {
        os << " " << _measurement[i];
    }

    for(int i = 0; i < _W.rows(); i++) {
        for(int j = 0; j < _W.cols(); j++) {
            os << " " << _W(i, j);
        }
    }

    for(int i = 0; i < _information.rows(); i++) {
        for(int j = i; j < _information.cols(); j++) {
            os << " " << _information(i, j);
        }
    }

    return true;
}
