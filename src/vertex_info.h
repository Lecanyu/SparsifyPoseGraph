/*
 * vertex_info.h
 *
 *  Created on: 14/oct/2014
 *      Author: Mladen Mazuran
 */

#ifndef VERTEX_INFO_H_
#define VERTEX_INFO_H_

#include <g2o/core/optimizable_graph.h>
#include <Eigen/Core>
#include <typeinfo>

class TypeInfo;

class VertexInfo {
public:
    virtual ~VertexInfo() {}

    virtual size_t id() const = 0;

    virtual std::string name() const = 0;

    virtual Eigen::VectorXd error(
            const g2o::HyperGraph::Vertex *vref,
            const g2o::HyperGraph::Vertex *vcmp) const = 0;

    virtual void setZero(g2o::HyperGraph::Vertex *v) const = 0;

    virtual void copyEstimate(
            const g2o::HyperGraph::Vertex *vsrc,
            g2o::HyperGraph::Vertex *vdst) const = 0;

    virtual g2o::HyperGraph::Vertex *create() const = 0;

    virtual bool isType(const g2o::HyperGraph::Vertex *v) const = 0;

    virtual g2o::HyperGraph::Vertex *clone(
            const g2o::HyperGraph::Vertex *v) const {
        g2o::HyperGraph::Vertex *vc = create();
        copyEstimate(v, vc);
        return vc;
    }

    virtual bool canReparametrize(const g2o::HyperGraph::Vertex *v) const = 0;

    virtual void reparametrize(const g2o::HyperGraph::Vertex *vsrc,
            g2o::HyperGraph::Vertex *vdst) const = 0;
};

template <typename T>
class VertexInfoBase : public VertexInfo
{
public:
    virtual ~VertexInfoBase() {}

    virtual size_t id() const {
        return typeid(T).hash_code();
    }

    virtual std::string name() const {
        return typeid(T).name();
    }

    virtual Eigen::VectorXd error(
            const g2o::HyperGraph::Vertex *vref,
            const g2o::HyperGraph::Vertex *vcmp) const {
        const T *vref2 = dynamic_cast<const T *>(vref);
        const T *vcmp2 = dynamic_cast<const T *>(vcmp);
        assert(vref2 && vcmp2 && "Wrong types");
        return error(vref2, vcmp2);
    }

    virtual Eigen::VectorXd error(const T *vref, const T *vcmp) const = 0;

    virtual void setZero(g2o::HyperGraph::Vertex *v) const {
        T *v2 = dynamic_cast<T *>(v);
        assert(v2 && "Wrong type");
        setZero(v2);
    }

    virtual void setZero(T *v) const = 0;

    virtual void copyEstimate(
            const g2o::HyperGraph::Vertex *vsrc,
            g2o::HyperGraph::Vertex *vdst) const {
        const T *vsrc2 = dynamic_cast<const T *>(vsrc);
        T *vdst2 = dynamic_cast<T *>(vdst);
        assert(vsrc2 && vdst2 && "Wrong types");
        copyEstimate(vsrc2, vdst2);
    }

    virtual void copyEstimate(const T *vsrc, T *vdst) const = 0;

    virtual g2o::HyperGraph::Vertex *create() const {
        return new T;
    }

    virtual bool isType(const g2o::HyperGraph::Vertex *v) const {
        return dynamic_cast<const T *>(v);
    }

    virtual bool canReparametrize(const g2o::HyperGraph::Vertex *v) const {
        const T *v2 = dynamic_cast<const T *>(v);
        return v2 && canReparametrize(v2);
    }

    virtual bool canReparametrize(const T *v) const = 0;

    virtual void reparametrize(const g2o::HyperGraph::Vertex *vsrc,
            g2o::HyperGraph::Vertex *vdst) const {
        const T *vsrc2 = dynamic_cast<const T *>(vsrc);
        T *vdst2 = dynamic_cast<T *>(vdst);
        assert(vsrc2 && vdst2 && "Wrong types");
        return reparametrize(vsrc2, vdst2);
    }

    virtual void reparametrize(const T *vsrc, T *vdst) const = 0;
};

class UnknownVertexInfo : public VertexInfo {
public:
    virtual size_t id() const { return typeid(UnknownVertexInfo).hash_code(); }

    virtual std::string name() const { return "Unknown"; }

    virtual Eigen::VectorXd error(
                const g2o::HyperGraph::Vertex *vref,
                const g2o::HyperGraph::Vertex *vcmp) const {
        const g2o::OptimizableGraph::Vertex *vref2 =
                dynamic_cast<const g2o::OptimizableGraph::Vertex *>(vref);
        if(vref2) {
            return Eigen::VectorXd::Zero(vref2->dimension());
        } else {
            return Eigen::VectorXd::Zero(0);
        }
    }

    virtual void setZero(g2o::HyperGraph::Vertex *) const {}

    virtual void copyEstimate(
            const g2o::HyperGraph::Vertex *,
            g2o::HyperGraph::Vertex *) const {}

    virtual g2o::HyperGraph::Vertex *create() const { return NULL; }

    virtual bool isType(const g2o::HyperGraph::Vertex *) const { return false; }

    virtual g2o::HyperGraph::Vertex *clone(
                const g2o::HyperGraph::Vertex *) const { return NULL; }

    virtual bool canReparametrize(const g2o::HyperGraph::Vertex *) const { return false; }

    virtual void reparametrize(const g2o::HyperGraph::Vertex *,
            g2o::HyperGraph::Vertex *) const {}
};

#endif /* VERTEX_INFO_H_ */
