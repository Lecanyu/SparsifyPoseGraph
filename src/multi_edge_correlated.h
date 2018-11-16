/*
 * multi_edge_correlated.h
 *
 *  Created on: 14/ott/2014
 *      Author: Mladen Mazuran
 */

#ifndef MULTI_EDGE_CORRELATED_H_
#define MULTI_EDGE_CORRELATED_H_

#include <g2o/core/base_multi_edge.h>
#include <Eigen/Core>
#include <list>

template <typename EdgeType>
class MultiEdgeCorrelated : public g2o::BaseMultiEdge<-1, std::list<typename EdgeType::Measurement> > {
public:
    typedef typename EdgeType::Measurement Measurement;
    typedef g2o::BaseMultiEdge<-1, std::list<typename EdgeType::Measurement> > Base;
    typedef std::list<std::list<int> > MappingsType;
    struct const_iterator;

    MultiEdgeCorrelated();
    virtual ~MultiEdgeCorrelated();

    void setDimension(int dimension_) {
        resize(0); //resize(dimension_);
        _dimension = dimension_;
        _information.resize(dimension_, dimension_);
        _error.resize(dimension_, 1);
        // _measurement.resize(dimension_, 1);
    }

    void setMeasurementCount(int measurements) {
        setDimension(measurements * EdgeType::Dimension);
    }

    int measurementCount() const {
        return _dimension / EdgeType::Dimension;
    }

    void addMeasurement(
            const g2o::HyperGraph::VertexContainer &verts,
            const Measurement &meas);


    virtual void computeError();

    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;

    virtual bool setMeasurementFromState();

    virtual void initialEstimate(
            const g2o::HyperGraph::VertexSet &fixed,
            g2o::OptimizableGraph::Vertex *toEstimate);
    virtual double initialEstimatePossible(
            const g2o::HyperGraph::VertexSet &fixed,
            g2o::OptimizableGraph::Vertex *toEstimate);

    virtual void linearizeOplus();

    const_iterator begin() const { return const_iterator(this); }
    const_iterator end() const { return const_iterator(this, true); }

    using Base::resize;

protected:
    using Base::_measurement;
    using Base::_information;
    using Base::_error;
    using Base::_vertices;
    using Base::_dimension;
    using Base::_jacobianOplus;

    g2o::HyperGraph::VertexContainer toVertexContainer(
            const std::list<int> &vmap) const;
    static g2o::HyperGraph::VertexSet setIntersection(
            const g2o::HyperGraph::VertexSet &a,
            const g2o::HyperGraph::VertexContainer &b);
    static bool setMembership(
            const g2o::HyperGraph::VertexContainer &a,
            const g2o::HyperGraph::Vertex *v);

public:
    struct const_iterator {
    public:
        const_iterator(const const_iterator &c) { *this = c; }
        ~const_iterator() {}

        const_iterator &operator=(const const_iterator &c) {
            _it = c._it; _itm = c._itm; _edge = c._edge; return *this;
        }

        const_iterator &operator++() {
           ++_it; ++_itm; return *this;
        }

        bool operator==(const const_iterator &c) const { return c._it == _it; }
        bool operator!=(const const_iterator &c) const { return c._it != _it; }
        const const_iterator &operator*() const { return *this; }

        g2o::HyperGraph::VertexContainer vertices() const {
            return _edge->toVertexContainer(*_it);
        }

        Measurement measurement() const {
            return *_itm;
        }

    private:
        const_iterator(const MultiEdgeCorrelated *edge, bool end = false) :
            _it(end ? edge->_mappings.end() : edge->_mappings.begin()),
            _itm(edge->_measurement.begin()), _edge(edge) {}

        MappingsType::const_iterator _it;
        typename std::list<Measurement>::const_iterator _itm;
        const MultiEdgeCorrelated *_edge;
        friend class MultiEdgeCorrelated;
    };

private:
    MappingsType _mappings;
    EdgeType *_mockedge;
    friend class const_iterator;
};


#include "multi_edge_correlated.hpp"

#endif /* MULTI_EDGE_CORRELATED_H_ */
