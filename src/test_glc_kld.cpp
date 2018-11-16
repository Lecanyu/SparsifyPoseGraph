/*
 * glc_kld_test.cpp
 *
 *  Created on: Oct 14, 2014
 *      Author: mazuran
 */

#include <isam/Slam.h>
#include <isam/slam3d.h>
#include <isam/glc.h>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <g2o/types/slam3d/isometry3d_mappings.h>
#include <g2o/stuff/misc.h>
#include <Eigen/Cholesky>

Eigen::MatrixXd readmat(std::istream &s, int r, int c)
{
    char skip;
    Eigen::MatrixXd m(r, c);
    for(int i = 0; i < r; i++) {
        for(int j = 0; j < c; j++) {
            s >> skip >> m(i,j);
        }
    }
    s >> skip;
    return m;
}

isam::Slam *load(const char *fname, std::map<uint64_t, isam::Pose3d_Node *> &lookup)
{
    isam::Slam *slam = new isam::Slam;
    isam::Properties prop;
    prop.verbose = true;
    prop.quiet = false;
    prop.mod_solve = isam::LEVENBERG_MARQUARDT;
    prop.epsilon_abs = 1e-8; //1e-5;
    prop.epsilon_rel = 1e-8; //1e-6;
    prop.mod_batch = 1;
    slam->set_properties(prop);

    std::ifstream f(fname);

    while(f.good()) {
        uint64_t id1, id2;
        char c;
        std::string line, tag;
        std::getline(f, line);
        std::stringstream ss(line);

        ss >> tag >> id1;

        if(tag == "NODE") {
            Eigen::VectorXd pose = readmat(ss, 6, 1);
            isam::Pose3d_Node *node = new isam::Pose3d_Node;
            node->init(isam::Pose3d(g2o::internal::fromVectorET(pose)));
            slam->add_node(node);
            lookup[id1] = node;
            // std::cout << pose(0) << " " << pose(1) << std::endl;
        } else if(tag == "FACTOR_PRIOR_POSE") {
            Eigen::VectorXd meas = readmat(ss, 6, 1);
            Eigen::MatrixXd cov = readmat(ss, 6, 6);
            isam::Pose3d_Factor *f = new isam::Pose3d_Factor(
                    lookup[id1], isam::Pose3d(g2o::internal::fromVectorET(meas)),
                    isam::Covariance(cov));
            cov.col(3).swap(cov.col(5));
            cov.row(3).swap(cov.row(5));
            slam->add_factor(f);
        } else if(tag == "FACTOR_POSE_POSE") {
            ss >> id2;
            Eigen::VectorXd meas = readmat(ss, 6, 1);
            Eigen::MatrixXd cov = readmat(ss, 6, 6);
            cov.col(3).swap(cov.col(5));
            cov.row(3).swap(cov.row(5));
            isam::Pose3d_Pose3d_Factor *f = new isam::Pose3d_Pose3d_Factor(
                    lookup[id1], lookup[id2], isam::Pose3d(g2o::internal::fromVectorET(meas)),
                    isam::Covariance(cov));
            slam->add_factor(f);
        }
    }

    f.close();

    slam->batch_optimization();
    slam->print_stats();

    return slam;
}

Eigen::VectorXd stack(isam::Slam *isam, std::list<isam::Node *> &keepList)
{
    Eigen::VectorXd ret(6 * keepList.size());
    int i = 0;
    for(isam::Node *n: keepList) {
        isam::Pose3d_Node *p3dn = dynamic_cast<isam::Pose3d_Node *>(n);
        ret.segment(i, p3dn->dim()) = p3dn->value().vector();
        i += p3dn->dim();
    }
    return ret;
}

int main(int argc, char *argv[])
{
    std::list<isam::Node *> keepList1, keepList2;
    std::map<uint64_t, isam::Pose3d_Node *> lookup1, lookup2;
    isam::Slam *slamOriginal = load(argv[1], lookup1);
    isam::Slam *slamMarginal = load(argv[1], lookup2);

    auto it1 = lookup1.begin(), it2 = lookup2.begin();
    for(int i = 0; i < lookup1.size(); i++) {
        if(i < 3 || (i + 1) % 4 == 0) {
            keepList1.push_back(it1->second);
            keepList2.push_back(it2->second);
        } else {
            isam::Node *toRemove = it2->second;
            std::vector<isam::Factor *> removed = glc_elim_factors(toRemove);
            glc_remove_node(*slamMarginal, toRemove, true, new isam::GLC_RootShift);

            delete toRemove;
            for(auto factor: removed) {
                delete factor;
            }
        }
        ++it1;
        ++it2;
    }

    slamMarginal->batch_optimization();

    Eigen::MatrixXd cov1 = slamOriginal->covariances().marginal(keepList1);
    Eigen::MatrixXd cov2 = slamMarginal->covariances().marginal(keepList2);
    Eigen::VectorXd m1 = stack(slamOriginal, keepList1);
    Eigen::VectorXd m2 = stack(slamMarginal, keepList2);
    Eigen::VectorXd d = m1 - m2;

    Eigen::MatrixXd info2 = cov2.ldlt().solve(Eigen::MatrixXd::Identity(d.rows(), d.rows()));
    info2 = 0.5 * (info2 + info2.transpose()).eval();

    for(int i = 0; i < d.rows(); i++) {
        if((i / 3) % 2 == 1) {
            d(i) = g2o::normalize_theta(d(i));
        }
    }

    double kld = 0.5 * (
            (info2 * cov1).trace() +
            d.transpose() * info2 * d - d.rows() -
            cov1.ldlt().vectorD().array().log().sum() +
            cov2.ldlt().vectorD().array().log().sum());

    std::cout << "KLD = " << kld << std::endl;

    return 0;
}

