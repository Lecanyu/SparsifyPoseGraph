/*
 * isometryxd.cpp
 *
 *  Created on: Jan 8, 2014
 *      Author: mazuran
 */

#include "isometryxd.h"
#include <g2o/types/slam3d/isometry3d_mappings.h>

IsometryXd::IsometryXd(bool is2d) :
        _se2(NULL), _se3(NULL)
{
    if(is2d) {
        _se2 = new g2o::SE2;
    } else {
        _se3 = new Eigen::Isometry3d;
        *_se3 = Eigen::Isometry3d::Identity();
    }
}

IsometryXd::IsometryXd(const g2o::SE2 &rt) :
        _se2(new g2o::SE2(rt)), _se3(NULL)
{
}

IsometryXd::IsometryXd(const isam::Pose2d &rt) :
        _se2(new g2o::SE2(rt.vector())), _se3(NULL)
{
}

IsometryXd::IsometryXd(const isam::Pose3d &rt) :
        _se2(NULL), _se3(new Eigen::Isometry3d(internal::fromVectorISAM(rt.vector())))
{
}

IsometryXd::IsometryXd(const Eigen::Isometry3d &rt) :
        _se2(NULL), _se3(new Eigen::Isometry3d(rt))
{
}

IsometryXd::IsometryXd(const Eigen::VectorXd &rt, IsometryXd::SE3Mode mode) :
        _se2(NULL), _se3(NULL)
{
    if(rt.rows() == 3) {
        _se2 = new g2o::SE2(rt);
    } else if(rt.rows() == 6) {
        if(mode == EulerAnglesISAM) {
            _se3 = new Eigen::Isometry3d(internal::fromVectorISAM(rt));
        } else if(mode == EulerAnglesG2O) {
            _se3 = new Eigen::Isometry3d(g2o::internal::fromVectorET(rt));
        } else if(mode == CondensedQuaternion) {
            _se3 = new Eigen::Isometry3d(g2o::internal::fromVectorMQT(rt));
        }
    } else /* if(rt.rows() == 7 && mode == Quaternion) */ {
        _se3 = new Eigen::Isometry3d(g2o::internal::fromVectorQT(rt));
    }
}

IsometryXd::IsometryXd(const IsometryXd &rt) :
        _se2(NULL), _se3(NULL)
{
    *this = rt;
}

IsometryXd::~IsometryXd()
{
    delete _se2;
    delete _se3;
}

IsometryXd IsometryXd::inverse() const
{
    if(_se2) {
        return IsometryXd(_se2->inverse());
    } else {
        return IsometryXd(_se3->inverse());
    }
}


IsometryXd &IsometryXd::operator=(const IsometryXd &rt)
{
    delete _se2;
    delete _se3;
    if(rt._se2) {
        _se2 = new g2o::SE2(*rt._se2);
        _se3 = NULL;
    } else {
        _se2 = NULL;
        _se3 = new Eigen::Isometry3d(*rt._se3);
    }
    return *this;
}

IsometryXd IsometryXd::operator*(const IsometryXd &rt) const
{
    if(_se2) {
        assert(rt._se2 && "Incompatible product");
        return IsometryXd((*_se2) * (*rt._se2));
    } else {
        assert(rt._se3 && "Incompatible product");
        return IsometryXd((*_se3) * (*rt._se3));
    }
}

IsometryXd &IsometryXd::operator*=(const IsometryXd &rt)
{
    return (*this) = (*this) * rt;
}

Eigen::VectorXd IsometryXd::vector(IsometryXd::SE3Mode mode) const
{
    if(_se2) {
        return _se2->toVector();
    } else {
        if(mode == EulerAnglesISAM) {
            return internal::toVectorISAM(*_se3);
        } else if(mode == EulerAnglesG2O) {
            return g2o::internal::toVectorET(*_se3);
        } else if(mode == CondensedQuaternion) {
            return g2o::internal::toVectorMQT(*_se3);
        } else /* if(mode == Quaternion) */ {
            return g2o::internal::toVectorQT(*_se3);
        }
    }
}

namespace internal {

Eigen::Matrix<double, 6, 1> toVectorISAM(const Eigen::Isometry3d &rt)
{
    Eigen::Matrix<double, 6, 1> ret = g2o::internal::toVectorET(rt);
    std::swap(ret[3], ret[5]);
    return ret;
}

Eigen::Isometry3d fromVectorISAM(const Eigen::Matrix<double, 6, 1> &rt)
{
    g2o::Vector6d rt2 = rt;
    std::swap(rt2[3], rt2[5]);
    return g2o::internal::fromVectorET(rt2);
}

} /* namespace internal */
