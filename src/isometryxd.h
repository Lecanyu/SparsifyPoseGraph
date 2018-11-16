/*
 * isometryxd.h
 *
 *  Created on: Jan 8, 2014
 *      Author: mazuran
 */

#ifndef ISOMETRYXD_H_
#define ISOMETRYXD_H_

#include <isam/Pose2d.h>
#include <isam/Pose3d.h>
#include <g2o/types/slam2d/se2.h>
#include <Eigen/Geometry>

class IsometryXd {
public:
    enum SE3Mode {
        EulerAnglesISAM, EulerAnglesG2O, CondensedQuaternion, Quaternion
    };

    IsometryXd(bool is2d = true);
    IsometryXd(const g2o::SE2 &rt);
    IsometryXd(const isam::Pose2d &rt);
    IsometryXd(const isam::Pose3d &rt);
    IsometryXd(const Eigen::Isometry3d &rt);
    IsometryXd(const Eigen::VectorXd &rt, SE3Mode mode = EulerAnglesISAM);
    IsometryXd(const IsometryXd &rt);
    virtual ~IsometryXd();

    IsometryXd inverse() const;

    IsometryXd &operator=(const IsometryXd &rt);
    IsometryXd operator*(const IsometryXd &rt) const;
    IsometryXd &operator*=(const IsometryXd &rt);

    Eigen::VectorXd vector(SE3Mode mode = EulerAnglesISAM) const;
    const g2o::SE2 &se2() const { return *_se2; }
    const Eigen::Isometry3d &se3() const { return *_se3; }

    bool is2d() const { return _se2; }
    bool is3d() const { return _se3; }

private:
    g2o::SE2 *_se2;
    Eigen::Isometry3d *_se3;
};

namespace internal {
    Eigen::Matrix<double, 6, 1> toVectorISAM(const Eigen::Isometry3d &rt);
    Eigen::Isometry3d fromVectorISAM(const Eigen::Matrix<double, 6, 1> &rt);
}

#endif /* ISOMETRYXD_H_ */
