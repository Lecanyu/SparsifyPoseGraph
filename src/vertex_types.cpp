/*
 * vertex_types.cpp
 *
 *  Created on: Oct 14, 2014
 *      Author: mazuran
 */

#include "type_info.h"
#include <g2o/types/slam2d/vertex_se2.h>
#include <g2o/types/slam3d/vertex_se3.h>
#include <g2o/types/slam3d/isometry3d_mappings.h>
#include "isometryxd.h"

class VertexSE2Info : public VertexInfoBase<g2o::VertexSE2> {
public:
    virtual Eigen::VectorXd error(
            const g2o::VertexSE2 *vref, const g2o::VertexSE2 *vcmp) const {
        Eigen::Vector3d d = vref->estimate().toVector() - vcmp->estimate().toVector();
        d(2) = g2o::normalize_theta(d(2));
        return d;
    }

    virtual void setZero(g2o::VertexSE2 *v) const {
        v->setEstimate(g2o::SE2(0, 0, 0));
    }

    virtual void copyEstimate(
            const g2o::VertexSE2 *vsrc, g2o::VertexSE2 *vdst) const {
        vdst->setEstimate(vsrc->estimate());
    }

    virtual bool canReparametrize(const g2o::VertexSE2 *) const {
        return true;
    }

    virtual void reparametrize(const g2o::VertexSE2 *vsrc, g2o::VertexSE2 *vdst) const {
        vdst->setEstimate(vsrc->estimate().inverse() * vdst->estimate());
    }

private:
    int _id;
    friend class TypeInfo;
};

class VertexSE3Info : public VertexInfoBase<g2o::VertexSE3> {
public:
    virtual Eigen::VectorXd error(
            const g2o::VertexSE3 *vref, const g2o::VertexSE3 *vcmp) const {
        // TODO: Which way?
        // Eigen::Isometry3d isod = vcmp->estimate().inverse() * vref->estimate();
        Eigen::Isometry3d isod = vref->estimate().inverse() * vcmp->estimate();
        Eigen::VectorXd d = g2o::internal::toVectorMQT(isod);
        return d;
    }

    virtual void setZero(g2o::VertexSE3 *v) const {
        v->setEstimate(Eigen::Isometry3d::Identity());
    }

    virtual void copyEstimate(
            const g2o::VertexSE3 *vsrc, g2o::VertexSE3 *vdst) const {
        vdst->setEstimate(vsrc->estimate());
    }

    virtual bool canReparametrize(const g2o::VertexSE3 *) const {
        return true;
    }

    virtual void reparametrize(const g2o::VertexSE3 *vsrc, g2o::VertexSE3 *vdst) const {
        vdst->setEstimate(vsrc->estimate().inverse() * vdst->estimate());
    }
private:
    int _id;
    friend class TypeInfo;
};

G2S_REGISTER_TYPE(VertexSE2Info);
G2S_REGISTER_TYPE(VertexSE3Info);

