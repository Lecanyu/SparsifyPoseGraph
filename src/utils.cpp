#include "utils.h"
#include <g2o/core/optimization_algorithm_with_hessian.h>
#include <isam/slam2d.h>
#include <isam/slam3d.h>
#include "isometryxd.h"
#include "se2_compatibility.h"
#include "se3_compatibility.h"

typedef VertexSE2ISAM VertexSE2Type;
typedef VertexSE3ISAM VertexSE3Type;

Eigen::VectorXd selectVariables(const Eigen::VectorXd &original, const std::vector<int> &variables)
{
    Eigen::VectorXd ret(variables.size());
    int i = 0;
    for(int vi: variables) {
        ret[i++] = original[vi];
    }
    return ret;
}

Eigen::MatrixXd selectVariables(const Eigen::MatrixXd &original, const std::vector<int> &variables)
{
    return selectVariables(original, variables, variables);
}

Eigen::MatrixXd selectVariables(const Eigen::MatrixXd &original, const std::vector<int> &variables1, const std::vector<int> &variables2)
{
    Eigen::MatrixXd ret(variables1.size(), variables2.size());
    int i = 0;
    for(int vi: variables1) {
        int j = 0;
        for(int vj: variables2) {
            ret(i, j) = original(vi, vj);
            j++;
        }
        i++;
    }
    return ret;
}

Eigen::SparseMatrix<double> selectVariables(const Eigen::SparseMatrix<double> &original,
        const std::vector<int> &variables)
{
    return selectVariables(original, variables, variables);
}

Eigen::SparseMatrix<double> selectVariables(const Eigen::SparseMatrix<double> &original,
        const std::vector<int> &variables1, const std::vector<int> &variables2)
{
    std::vector<Eigen::Triplet<double> > trs;
    for(int k = 0; k < original.outerSize(); ++k) {
        for(Eigen::SparseMatrix<double>::InnerIterator it(original, k); it; ++it) {
            auto x = std::lower_bound(variables1.begin(), variables1.end(), it.row());
            if(x != variables1.end() && *x == it.row()) {
                auto y = std::lower_bound(variables2.begin(), variables2.end(), it.col());
                if(y != variables2.end() && *y == it.col()) {
                    trs.push_back(Eigen::Triplet<double>(
                            x - variables1.begin(), y - variables2.begin(), it.value()));
                }
            }
        }
    }

    Eigen::SparseMatrix<double> ret(variables1.size(), variables2.size());
    ret.setFromTriplets(trs.begin(), trs.end());
    return ret;
}

double kullbackLeiblerDivergence(
        const Eigen::VectorXd &diff, const Eigen::MatrixXd &infox,
        const Eigen::MatrixXd &maty, KullbackMode mode)
{
    Eigen::LDLT<Eigen::MatrixXd> ldlty(maty);
    double logdetx = infox.ldlt().vectorD().array().log().sum(), logdety, innerprod;
    if(mode == InformationCovariance) {
        logdety = ldlty.vectorD().array().log().sum();
        innerprod = infox.cwiseProduct(maty).sum();
    } else {
        logdety = - ldlty.vectorD().array().log().sum();
        innerprod = ldlty.solve(infox).trace();
    }

    double mahalanobis = diff.transpose() * infox * diff;
    double kld = 0.5 * (innerprod + mahalanobis - logdetx - logdety - infox.rows());

//    if(kld > 1000) {
//        std::cout << "diff        = [" << diff.transpose() << "];" << std::endl;
//        std::cout << "logdetx     = " << logdetx << std::endl;
//        std::cout << "logdety     = " << logdety << std::endl;
//        std::cout << "innerprod   = " << innerprod << std::endl;
//        std::cout << "mahalanobis = " << mahalanobis << std::endl;
//        std::cout << "squarednorm = " << diff.squaredNorm() << std::endl << std::endl;
//        //if(kld > 0.4) exit(0);
//    }
    return kld;
}

std::vector<Eigen::Triplet<double> > triplets(const g2o::SparseBlockMatrix<Eigen::MatrixXd> &mat)
{
    std::vector<Eigen::Triplet<double> > ret;
    const std::vector<g2o::SparseBlockMatrix<Eigen::MatrixXd>::IntBlockMap> &blockCols = mat.blockCols();
    for(size_t c = 0; c < blockCols.size(); ++c) {
        for(g2o::SparseBlockMatrix<Eigen::MatrixXd>::IntBlockMap::const_iterator it = blockCols[c].begin();
                it != blockCols[c].end(); ++it) {

            const int r = it->first;
            const Eigen::MatrixXd &m = *(it->second);

            for(int cc = 0; cc < m.cols(); ++cc) {
                for(int rr = 0; rr < m.rows(); ++rr) {
                    int aux_r = mat.rowBaseOfBlock(r) + rr;
                    int aux_c = mat.colBaseOfBlock(c) + cc;
                    ret.push_back(Eigen::Triplet<double>(aux_r, aux_c, m(rr, cc)));
                    if(r != c) {
                        ret.push_back(Eigen::Triplet<double>(aux_c, aux_r, m(rr, cc)));
                    }
                }
            }
        }
    }
    return ret;
}

std::vector<Eigen::Triplet<double> > triplets(const Eigen::SparseMatrix<double> &mat)
{
    std::vector<Eigen::Triplet<double> > ret;
    for (int k = 0; k < mat.outerSize(); k++) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
            ret.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
        }
    }
    return ret;
}

Eigen::SparseMatrix<double> sparse(const isam::SparseSystem &mat)
{
    Eigen::SparseMatrix<double> ret(mat.num_rows(), mat.num_cols());
    const int *lookup = mat.r_to_a();
    std::vector<Eigen::Triplet<double> > trs;
    for(int r = 0; r < mat.num_rows(); r++) {
        const isam::SparseVector &row = mat.get_row(r);
        isam::SparseVectorIter svi(row);
        while(svi.valid()) {
            trs.push_back(Eigen::Triplet<double>(r, lookup[svi.get()], svi.get_val()));
            svi.next();
        }
    }
    ret.setFromTriplets(trs.begin(), trs.end());
    return ret;
}

Eigen::SparseMatrix<double> sparse(const g2o::SparseBlockMatrix<Eigen::MatrixXd> &mat)
{
    Eigen::SparseMatrix<double> ret(mat.rows(), mat.cols());
    std::vector<Eigen::Triplet<double> > trs = triplets(mat);
    ret.setFromTriplets(trs.begin(), trs.end());
    return ret;
}

Eigen::MatrixXd dense(const g2o::SparseBlockMatrix<Eigen::MatrixXd> &mat)
{
    Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(mat.rows(), mat.cols());
    std::vector<Eigen::Triplet<double> > trs = triplets(mat);
    for(auto &t: trs) {
        ret(t.row(), t.col()) = t.value();
    }
    return ret;
}

int fullSize(g2o::SparseOptimizer *so)
{
    int fullsize = 0;
    /**
     * Canyu Le
     * We directly ignore one vertex instead of hardcoding vertex Id=0
     * */
    // -----------------------
    bool ignore_start_v = true;
    for(auto v: so->vertices()) {
        if(ignore_start_v)
            ignore_start_v = false;
        else
            fullsize += static_cast<g2o::OptimizableGraph::Vertex *>(v.second)->dimension();
    }
    // -----------------------


    return fullsize;
}

int fullSize(isam::Slam *isam)
{
    int fullsize = 0;
    for(auto n: isam->get_nodes()) {
        fullsize += n->dim();
    }
    return fullsize;
}

Eigen::VectorXd stack(g2o::SparseOptimizer *so)
{
    Eigen::VectorXd ret(fullSize(so));
    int i = 0;
    for(auto v: so->indexMapping()) {
        VertexSE2Type *vse2 = dynamic_cast<VertexSE2Type *>(v);
        VertexSE3Type *vse3 = dynamic_cast<VertexSE3Type *>(v);
        assert((vse2 || vse3) && "Can handle only SE2 or SE3 vertices for now");
        if(vse2) {
            ret.block(i, 0, v->dimension(), 1) =
                    vse2->estimate().toVector();
        } else if(vse3) {
            ret.block(i, 0, v->dimension(), 1) =
                    //g2o::internal::toVectorMQT(vse3->estimate());
                    internal::toVectorISAM(vse3->estimate());
        }
        i += v->dimension();
    }
    return ret;
}

Eigen::VectorXd stack(isam::Slam *isam)
{
    Eigen::VectorXd ret(fullSize(isam));
    int i = 0;
    for(auto n: isam->get_nodes()) {
        isam::Pose2d_Node *p2dn = dynamic_cast<isam::Pose2d_Node *>(n);
        isam::Pose3d_Node *p3dn = dynamic_cast<isam::Pose3d_Node *>(n);
        assert((p2dn || p3dn) && "Can handle only Pose2d or Pose3d nodes for now");
        if(p2dn) {
            ret.block(i, 0, p2dn->dim(), 1) = p2dn->value().vector();
            i += p2dn->dim();
        } else if(p3dn) {
            ret.block(i, 0, p3dn->dim(), 1) = p3dn->value().vector();
            i += p3dn->dim();
        }
    }
    return ret;
}
