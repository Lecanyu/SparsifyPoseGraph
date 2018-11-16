/*
 * se3_compatibility.h
 *
 *  Created on: 11/gen/2014
 *      Author: Mladen Mazuran
 */

#ifndef SE3_COMPATIBILITY_H_
#define SE3_COMPATIBILITY_H_

#include <g2o/types/slam3d/edge_se3.h>
#include <g2o/types/slam3d/isometry3d_mappings.h>
#include <g2o/types/slam3d/isometry3d_gradients.h>
#include <Eigen/Geometry>
#include "isometryxd.h"

typedef g2o::VertexSE3 VertexSE3ISAM;

/*
 * The implementation of g2o::EdgeSE3 works with condensed quaternions with the difference being
 * computed on the manifold. On the other hand, isam works with Euler angles and the difference
 * is computed in the Euclidean space (with angle wrapping). EdgeSE3ISAM reimplements this
 * behavior in order to provide a fair comparison between g2o and isam.
 */
class EdgeSE3ISAM : public g2o::EdgeSE3
{
public:
    EdgeSE3ISAM() : g2o::EdgeSE3() {}

#ifndef G2S_QUATERNIONS
    void computeError() {
        VertexSE3ISAM *from = static_cast<VertexSE3ISAM *>(_vertices[0]);
        VertexSE3ISAM *to   = static_cast<VertexSE3ISAM *>(_vertices[1]);
        _error = internal::toVectorISAM(from->estimate().inverse() * to->estimate()) -
                internal::toVectorISAM(_measurement);
        _error[3] = g2o::normalize_theta(_error[3]);
        _error[4] = g2o::normalize_theta(_error[4]);
        _error[5] = g2o::normalize_theta(_error[5]);
    }

#define G2S_SE3_CLOSED_FORM
#ifdef G2S_SE3_CLOSED_FORM
    void linearizeOplus() {
        VertexSE3ISAM *from = static_cast<VertexSE3ISAM *>(_vertices[0]);
        VertexSE3ISAM *to   = static_cast<VertexSE3ISAM *>(_vertices[1]);

        /*
         * g2o::EdgeSE3 computes the gradient of z^-1 * xi^-1 * xj in condensed quaternion form
         * wrt condensed quaternions. We need xi^-1 * xj - z in Euler angles.
         * Rather than implementing it from scratch, we can set z = identity, compute the gradient
         * with g2o and then map the condensed quaternions to Euler angles by using the chain
         * rule. This is probably less efficient than directly computing the gradient, but avoids
         * writing a boatload of formulas, which is fine by me.
         */
        Eigen::Isometry3d E;
        const Eigen::Isometry3d& Xi = from->estimate();
        const Eigen::Isometry3d& Xj = to->estimate();
        const Eigen::Isometry3d  Z  = Eigen::Isometry3d::Identity();

        /* With this we get the jacobian as computed by g2o */
        g2o::internal::computeEdgeSE3Gradient(E, _jacobianOplusXi , _jacobianOplusXj, Z, Xi, Xj);

        /* Map the measure in condensed quaternion form to the measure in Euler angles form */
        g2o::Vector6d eulerOminus = internal::toVectorISAM(E);
        Eigen::Matrix3d mapQuatToEulerOminus =
                jacobianQuatWRTEuler(eulerOminus.tail<3>()).inverse();
        Eigen::Matrix<double, 3, 6> JXi = mapQuatToEulerOminus * _jacobianOplusXi.bottomRows<3>();
        Eigen::Matrix<double, 3, 6> JXj = mapQuatToEulerOminus * _jacobianOplusXj.bottomRows<3>();
        _jacobianOplusXi.bottomRows<3>() = JXi;
        _jacobianOplusXj.bottomRows<3>() = JXj;
    }
#else /* G2S_SE3_CLOSED_FORM */
    void linearizeOplus() {
        std::cout << "diocan" << std::endl;
        const double eps = 0.0001;
        VertexSE3ISAM *verts[] = {
                static_cast<VertexSE3ISAM *>(_vertices[0]),
                static_cast<VertexSE3ISAM *>(_vertices[1])
        };

        g2o::Vector6d errorBefore = _error;
        double delta[6];
        std::fill(delta, delta + 6, 0.0);
        for(int i = 0; i < 2; i++) {
            VertexSE3ISAM *vert = verts[i];
            for(int j = 0; j < 6; j++) {
                vert->push();
                delta[j] = eps;
                vert->oplus(delta);
                computeError();
                g2o::Vector6d err = _error;
                vert->pop();
                
                vert->push();
                delta[j] = -eps;
                vert->oplus(delta);
                computeError();
                err -= _error;
                vert->pop();
                
                delta[j] = 0;
                if(i == 0) {
                    _jacobianOplusXi.col(j) = err / eps;
                } else {
                    _jacobianOplusXj.col(j) = err / eps;
                }
            }
        }
        _error = errorBefore;      
    }

#endif /* G2S_SE3_CLOSED_FORM */

#endif /* G2S_QUATERNIONS */

    /* Jacobian of condensed quaternion wrt Euler angles in isam order (yaw,pitch,roll) */
    template <typename Derived>
    static Eigen::Matrix3d jacobianQuatWRTEuler(const Eigen::MatrixBase<Derived> &v) {
        Eigen::Matrix3d ret;
        const double y2 = 0.5 * v[0], p2 = 0.5 * v[1], r2 = 0.5 * v[2];
        const double cy = std::cos(y2), sy = std::sin(y2);
        const double cp = std::cos(p2), sp = std::sin(p2);
        const double cr = std::cos(r2), sr = std::sin(r2);

        ret <<
                0.5 * (-sr * cp * sy - cr * sp * cy),   /* dq1/dy */
                0.5 * (-sr * sp * cy - cr * cp * sy),   /* dq1/dp */
                0.5 * ( cr * cp * cy + sr * sp * sy),   /* dq1/dr */

                0.5 * (-cr * sp * sy + sr * cp * cy),   /* dq2/dy */
                0.5 * ( cr * cp * cy - sr * sp * sy),   /* dq2/dp */
                0.5 * (-sr * sp * cy + cr * cp * sy),   /* dq2/dr */

                0.5 * ( cr * cp * cy + sr * sp * sy),   /* dq3/dy */
                0.5 * (-cr * sp * sy - sr * cp * cy),   /* dq3/dp */
                0.5 * (-sr * cp * sy - cr * sp * cy);   /* dq3/dr */
        return ret;
    }

    /*
     * Convert an information matrix expressed in condensed quaternion form to one expressed in
     * Euler angles form.
     */
    template <typename Derived>
    static Eigen::Matrix<double, 6, 6> eulerInformationFromQuat(
            const Eigen::MatrixBase<Derived> &info) {
        /*
         * jacobianQuatWRTEuler evaluated at zero is an anti-diagonal matrix of value 1/2.
         * To propagate the covariance we need the inverse of this, and putting it in information
         * form we invert it again. This means that to get the final result we have to multiply
         * the input information left and right by 0.5*I and then swap yaw with roll.
         */
        Eigen::Matrix<double, 6, 6> ret;
        ret.topLeftCorner<3, 3>()     = info.template topLeftCorner<3, 3>();
        ret.topRightCorner<3, 3>()    = 0.5 * info.template topRightCorner<3, 3>();
        ret.bottomLeftCorner<3, 3>()  = 0.5 * info.template bottomLeftCorner<3, 3>();
        ret.bottomRightCorner<3, 3>() = 0.25 * info.template bottomRightCorner<3, 3>();

        ret.col(3).swap(ret.col(5));
        ret.row(3).swap(ret.row(5));

        return ret;
    }

    /*
     * Convert an information matrix expressed in Euler angles form to one expressed in
     * condensed quaternions.
     */
    template <typename Derived>
    static Eigen::Matrix<double, 6, 6> quatInformationFromEuler(
            const Eigen::MatrixBase<Derived> &info) {
        Eigen::Matrix<double, 6, 6> ret;
        ret.topLeftCorner<3, 3>()     = info.template topLeftCorner<3, 3>();
        ret.topRightCorner<3, 3>()    = 2 * info.template topRightCorner<3, 3>();
        ret.bottomLeftCorner<3, 3>()  = 2 * info.template bottomLeftCorner<3, 3>();
        ret.bottomRightCorner<3, 3>() = 4 * info.template bottomRightCorner<3, 3>();

        ret.col(3).swap(ret.col(5));
        ret.row(3).swap(ret.row(5));

        return ret;
    }
};

#endif /* SE3_COMPATIBILITY_H_ */
