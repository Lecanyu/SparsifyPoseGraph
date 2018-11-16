/*
 * edge_info.h
 *
 *  Created on: 14/oct/2014
 *      Author: Mladen Mazuran
 */

#ifndef EDGE_INFO_H_
#define EDGE_INFO_H_

#include <g2o/core/optimizable_graph.h>
#include <typeinfo>

class TypeInfo;

class EdgeInfo {
public:
    virtual ~EdgeInfo() {}

    virtual size_t id() const = 0;

    virtual std::string name() const = 0;

    virtual g2o::HyperGraph::Edge *create() const = 0;

    virtual bool isType(const g2o::HyperGraph::Edge *e) const = 0;

    virtual g2o::HyperGraph::Edge *clone(const g2o::HyperGraph::Edge *e) const = 0;
};

template <typename T>
class EdgeInfoBase : public EdgeInfo
{
public:
    virtual ~EdgeInfoBase() {}

    virtual size_t id() const {
        return typeid(T).hash_code();
    }

    virtual std::string name() const {
        return typeid(T).name();
    }

    virtual bool isType(const g2o::HyperGraph::Edge *e) const {
        return dynamic_cast<const T *>(e);
    }

    virtual g2o::HyperGraph::Edge *create() const {
        return new T;
    }

    virtual g2o::HyperGraph::Edge *clone(const g2o::HyperGraph::Edge *e) const {
        const T *e2 = dynamic_cast<const T *>(e);
        assert(e2 && "Wrong type");
        return clone(e2);
    }

    virtual g2o::HyperGraph::Edge *clone(const T *e) const {
        T *c = new T;

        c->vertices().resize(e->vertices().size());

        for(int i = 0; i < e->vertices().size(); i++) {
            c->setVertex(i, e->vertices()[i]);
        }

        c->setMeasurement(e->measurement());
        c->setInformation(e->information());
        return c;
    }
};

class UnknownEdgeInfo : public EdgeInfo {
public:
    virtual size_t id() const { return typeid(UnknownEdgeInfo).hash_code(); }

    virtual std::string name() const { return "Unknown"; }

    virtual g2o::HyperGraph::Edge *create() const { return NULL; }

    virtual bool isType(const g2o::HyperGraph::Edge *e) const { return false; }

    virtual g2o::HyperGraph::Edge *clone(const g2o::HyperGraph::Edge *e) const { return NULL; }
};




#endif /* EDGE_INFO_H_ */
