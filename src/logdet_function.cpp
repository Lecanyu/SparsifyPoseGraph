/*
 * logdetfunction.cpp
 *
 *  Created on: May 22, 2014
 *      Author: mazuran
 */

#include "logdet_function.h"
#include <Eigen/Eigen>
#include <iostream>
#include <cfloat>
#include <cmath>

LogdetFunction::LogdetFunction(
        const JacobianMapping &mapping,
        const Eigen::MatrixXd &targetInformation) :
        _mapping(mapping), _target(targetInformation), _jacsize(0)
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(_target);

    for(auto &jaclist: _mapping) {
        _jacsize += jaclist.front().first.rows();
    }

#ifdef G2S_MORE_GENERIC_LESS_WORKING
    double maxeig = eig.eigenvalues().maxCoeff();
    double eps = std::nextafter(maxeig, DBL_MAX) - maxeig;
    double tol = std::max(_target.cols() * eps, 1e-8);

    int zeroeigs = (eig.eigenvalues().array() < tol).count();

    _S = eig.eigenvalues().tail(_target.rows() - zeroeigs).cwiseInverse();
    _U = eig.eigenvectors().rightCols(_target.rows() - zeroeigs);
#else /* G2S_MORE_GENERIC_LESS_WORKING */
    static const double cutoff = 1e-5;
    int smalleigs = (eig.eigenvalues().array() < cutoff).count();
    int dim = _mapping.front().front().first.cols();

    if(smalleigs <= dim) {
        _S = eig.eigenvalues().tail(_target.rows() - dim).cwiseInverse();
        _U = eig.eigenvectors().rightCols(_target.rows() - dim);
    } else {
        Eigen::MatrixXd candidates(_target.rows(), smalleigs);

        for(int i = 0; i < candidates.cols(); i++) {
            candidates.col(i) = eig.eigenvectors().col(i);
        }

        std::set<int> toDrop = chooseDimensions(candidates);
        _S = Eigen::VectorXd(_target.rows() - dim);
        _U = Eigen::MatrixXd(_target.rows(), _target.rows() - dim);

        for(int i = 0, j = 0; i < _target.rows(); i++) {
            if(toDrop.count(i) == 0) {
                _S[j] = std::min(std::abs(1 / eig.eigenvalues()[i]), 1e6 / eig.eigenvalues()[_target.rows() - 1]);
                _U.col(j) = eig.eigenvectors().col(i);
                j++;
            }
        }
    }
#endif /* G2S_MORE_GENERIC_LESS_WORKING */

    _logdet = _S.array().log().sum();
}

std::set<int> LogdetFunction::chooseDimensions(const Eigen::MatrixXd &candidates) const
{
    int dim = _mapping.front().front().first.cols();

    Eigen::MatrixXd JU = sparseJacobian() * candidates;
    std::vector<std::pair<double, int> > norms;
    for(int i = 0; i < JU.cols(); i++) {
        norms.push_back(std::make_pair(JU.col(i).norm(), i));
    }
    std::sort(norms.begin(), norms.end());
    std::set<int> toDrop;
    for(int i = 0; i < dim; i++) {
        toDrop.insert(norms[i].second);
    }
    return toDrop;
}

bool LogdetFunction::hasClosedFormSolution() const
{
    return _jacsize == _S.rows();
}

std::list<Eigen::MatrixXd> LogdetFunction::decondense(const Eigen::VectorXd &x) const
{
    std::list<Eigen::MatrixXd> ret;
    int k = 0;
    for(auto &jaclist: _mapping) {
        int n = jaclist.front().first.rows();
        Eigen::Map<const Eigen::MatrixXd> m(x.data() + k, n, n);
        k += n * n;
        ret.push_back(m.selfadjointView<Eigen::Lower>());
    }
    return ret;
}

Eigen::VectorXd LogdetFunction::condense(const std::list<Eigen::MatrixXd> &X) const
{
    int vecsize = 0, k = 0;
    for(const Eigen::MatrixXd &entry: X) {
        vecsize += entry.rows() * entry.cols();
    }

    Eigen::VectorXd ret(vecsize);

    for(const Eigen::MatrixXd &entry: X) {
        int entrysize = entry.rows() * entry.cols();
        Eigen::Map<const Eigen::VectorXd> m(entry.data(), entrysize);
        ret.segment(k, entrysize) = m;
        k += entrysize;
    }
    return ret;
}

double LogdetFunction::value(const Eigen::VectorXd &x)
{
    Eigen::MatrixXd UJXJU = _U.transpose() * informationProduct(x) * _U;
    UJXJU.triangularView<Eigen::StrictlyLower>() = UJXJU.triangularView<Eigen::StrictlyUpper>().transpose();
    _chol = Eigen::LDLT<Eigen::MatrixXd>(UJXJU);
    if(_chol.isPositive() && (_chol.vectorD().array() > 0).all()) {
        double kld = 0.5 * (
                (UJXJU.diagonal().array() * _S.array()).sum() -
                _chol.vectorD().array().log().sum() - _logdet - _S.rows());

        return kld;
    } else {
        return INFINITY;
    }
}

void LogdetFunction::gradient(const Eigen::VectorXd &x, Eigen::VectorXd &g)
{
    if(_chol.isPositive() && (_chol.vectorD().array() > 0).all()) {
        _xinv = _chol.solve(Eigen::MatrixXd::Identity(_S.rows(), _S.rows()));
        Eigen::MatrixXd middle = -_xinv;
        middle += _S.asDiagonal();
        Eigen::MatrixXd Y = _U * middle * _U.transpose();

        int k = 0;
        g.resize(x.rows());
        for(auto &jaclist: _mapping) {
            int n = jaclist.front().first.rows();
            Eigen::MatrixXd block = Eigen::MatrixXd::Zero(n, n);
            for(MeasurementJacobian::const_iterator sweep1 = jaclist.begin();
                    sweep1 != jaclist.end(); ++sweep1) {
                int m = sweep1->first.cols();
                Eigen::MatrixXd thisBlock =
                        sweep1->first *
                        Y.block(sweep1->second, sweep1->second, m, m) *
                        sweep1->first.transpose();
                block += 0.5 * (thisBlock + thisBlock.transpose());
                MeasurementJacobian::const_iterator sweep2 = sweep1;
                for(++sweep2; sweep2 != jaclist.end(); ++sweep2) {
                    int p = sweep2->first.cols();
                    thisBlock =
                            sweep1->first *
                            Y.block(sweep1->second, sweep2->second, m, p) *
                            sweep2->first.transpose();
                    block += thisBlock + thisBlock.transpose();
                }

            }

            Eigen::Map<Eigen::VectorXd> gsegment(block.data(), n * n);
            g.segment(k, n * n) = 0.5 * gsegment;
            k += n * n;
        }
    } else {
        g = Eigen::VectorXd::Zero(x.rows());
    }
}

void LogdetFunction::hessian(const Eigen::VectorXd &x, Eigen::MatrixXd &H)
{
    Eigen::MatrixXd JU = sparseJacobian() * _U;
    Eigen::MatrixXd P = JU * _xinv * JU.transpose();
    P = (0.5 * (P + P.transpose())).eval();

    int fullsize = 0, s = 0, k = 0;
    for(auto &jaclist: _mapping) {
        int n = jaclist.front().first.rows();
        fullsize += n * n;
    }

    Eigen::VectorXd hij(fullsize);
    for(auto &jaclist: _mapping) {
        int n = jaclist.front().first.rows();
        for(int jj = 0; jj < n; jj++) {
            for(int ii = 0; ii < n; ii++, s++) {
                singleVariableHessian(P, k + ii, k + jj, hij);
                H.row(s) = hij.transpose();
            }
        }
        k += n;
    }
}

void LogdetFunction::singleVariableHessian(const Eigen::MatrixXd &P, int i, int j, Eigen::VectorXd &hij) const
{
    int k = 0, s = 0;
    for(auto &jaclist: _mapping) {
        int n = jaclist.front().first.rows();
        for(int vv = 0; vv < n; vv++) {
            for(int uu = 0; uu < n; uu++) {
                hij(s++) = P(k + uu, i) * P(j, k + vv);
            }
        }
        k += n;
    }
}

Eigen::VectorXd LogdetFunction::educatedGuess()
{
    int fullsize = 0, k = 0;
    for(auto &jaclist: _mapping) {
        int n = jaclist.front().first.rows();
        fullsize += n * n;
    }

    Eigen::VectorXd ret = Eigen::VectorXd::Zero(fullsize);

    for(auto &jaclist: _mapping) {
        int n = jaclist.front().first.rows();
        Eigen::Map<Eigen::MatrixXd> X(ret.data() + k, n, n);
        X = Eigen::MatrixXd::Identity(n, n);
        k += n * n;
    }

    return ret;
}

std::list<Eigen::MatrixXd> LogdetFunction::closedFormSolution() const
{
    std::list<Eigen::MatrixXd> ret;
    Eigen::MatrixXd Sigma = _U * _S.asDiagonal() * _U.transpose();
    Sigma.triangularView<Eigen::StrictlyUpper>() = Sigma.triangularView<Eigen::StrictlyLower>().transpose();

    /* Check if we have a dense solution, if so use sparse jacobian, it's more efficient */
    if(_mapping.size() == 1) {
        Eigen::SparseMatrix<double> J = sparseJacobian();
        Eigen::MatrixXd piece = Sigma * J.transpose();
        Eigen::LLT<Eigen::MatrixXd> chol(J * piece);
        ret.push_back(chol.solve(Eigen::MatrixXd::Identity(J.rows(), J.rows())));
    } else {
        for(auto &jaclist: _mapping) {
            int n = jaclist.front().first.rows();
            Eigen::MatrixXd block = Eigen::MatrixXd::Zero(n, n);
            for(MeasurementJacobian::const_iterator sweep1 = jaclist.begin();
                    sweep1 != jaclist.end(); ++sweep1) {
                int m = sweep1->first.cols();

                Eigen::MatrixXd thisBlock =
                        sweep1->first *
                        Sigma.block(sweep1->second, sweep1->second, m, m) *
                        sweep1->first.transpose();
                block += 0.5 * (thisBlock + thisBlock.transpose());
                MeasurementJacobian::const_iterator sweep2 = sweep1;
                for(++sweep2; sweep2 != jaclist.end(); ++sweep2) {
                    int k = sweep2->first.cols();
                    thisBlock =
                            sweep1->first *
                            Sigma.block(sweep1->second, sweep2->second, m, k) *
                            sweep2->first.transpose();
                    block += thisBlock + thisBlock.transpose();
                }
            }

            /* TODO: Check if Cholesky doesn't have numerical issues for badly conditioned matrices */
            Eigen::LLT<Eigen::MatrixXd> chol(block);
            ret.push_back(chol.solve(Eigen::MatrixXd::Identity(n, n)));
        }
    }
    return ret;

}

Eigen::MatrixXd LogdetFunction::informationProduct(const Eigen::VectorXd &x) const
{
    std::list<Eigen::MatrixXd> X = decondense(x);
    return informationProduct(X);
}

Eigen::MatrixXd LogdetFunction::informationProduct(const std::list<Eigen::MatrixXd> &X) const
{
    Eigen::MatrixXd JXJ = Eigen::MatrixXd::Zero(_target.cols(), _target.rows());

    /* TODO: The sparse product might actually be always more efficient, need to check this */
    /* Check if we have a dense solution, if so use sparse jacobian, it's more efficient */
    if(X.size() == 1) {
        Eigen::SparseMatrix<double> J = sparseJacobian();
        Eigen::MatrixXd piece = J.transpose() * (*X.begin());
        JXJ = piece * J;
    } else {
        int k = 0;
        std::list<Eigen::MatrixXd>::const_iterator it = X.begin();
        for(auto &jaclist: _mapping) {
            int n = jaclist.front().first.rows();
            for(MeasurementJacobian::const_iterator sweep1 = jaclist.begin();
                    sweep1 != jaclist.end(); ++sweep1) {
                for(MeasurementJacobian::const_iterator sweep2 = sweep1;
                        sweep2 != jaclist.end(); ++sweep2) {
                    if(sweep1->second <= sweep2->second) {
                        JXJ.block(
                                sweep1->second, sweep2->second,
                                sweep1->first.cols(), sweep2->first.cols()) +=
                                sweep1->first.transpose() * (*it) * sweep2->first;
                    } else {
                        JXJ.block(
                                sweep2->second, sweep1->second,
                                sweep2->first.cols(), sweep1->first.cols()) +=
                                sweep2->first.transpose() * (*it) * sweep1->first;
                    }
                }
            }
            ++it;
        }
    }
    return JXJ.selfadjointView<Eigen::Upper>();
}

Eigen::SparseMatrix<double> LogdetFunction::sparseJacobian() const
{
    Eigen::SparseMatrix<double> J(_jacsize, _target.cols());
    std::vector<Eigen::Triplet<double> > triplets;
    int k = 0;
    for(auto &jaclist: _mapping) {
        for(auto &Jpair: jaclist) {
            for(int ii = 0; ii < Jpair.first.rows(); ii++) {
                for(int jj = 0; jj < Jpair.first.cols(); jj++) {
                    /* TODO: epsilon() is about 2.22045e-16, check if this is too big */
                    if(std::abs(Jpair.first(ii, jj)) >= std::numeric_limits<double>::epsilon()) {
                        triplets.push_back(Eigen::Triplet<double>(
                                k + ii, Jpair.second + jj, Jpair.first(ii, jj)));
                    }
                }
            }
        }
        k += jaclist.front().first.rows();
    }
    J.setFromTriplets(triplets.begin(), triplets.end());
    return J;
}

LogdetFunctionWithConstraints::LogdetFunctionWithConstraints(
        const JacobianMapping &mapping,
        const Eigen::MatrixXd &targetInformation) :
        LogdetFunction(mapping, targetInformation), _rho(0)
{
}

double LogdetFunctionWithConstraints::value(const Eigen::VectorXd &x)
{
    double original = LogdetFunction::value(x);
    std::list<Eigen::MatrixXd> X = decondense(x);
    Eigen::LDLT<Eigen::MatrixXd> chol;
    for(const Eigen::MatrixXd &Xblock: X) {
        chol.compute(Xblock);
        if(chol.isPositive() && (chol.vectorD().array() > 0).all()) {
            original -= _rho * chol.vectorD().array().log().sum();
        } else {
            return INFINITY;
        }
    }

    return original;
}

void LogdetFunctionWithConstraints::gradient(const Eigen::VectorXd &x, Eigen::VectorXd &g)
{
    LogdetFunction::gradient(x, g);

    std::list<Eigen::MatrixXd> X = decondense(x);
    Eigen::LDLT<Eigen::MatrixXd> chol;
    int k = 0;
    _invXblocks.clear();
    if(_chol.isPositive() && (_chol.vectorD().array() > 0).all()) {
        for(const Eigen::MatrixXd &Xblock: X) {
            chol.compute(Xblock);
            if(chol.isPositive() && (chol.vectorD().array() > 0).all()) {
                int n = Xblock.rows();
                Eigen::MatrixXd invXblock = chol.solve(Eigen::MatrixXd::Identity(n, n));
                _invXblocks.push_back(invXblock);
                Eigen::Map<Eigen::VectorXd> gsegment(invXblock.data(), n * n);
                g.segment(k, n * n) -= _rho * gsegment;
                k += n * n;
            } else {
                g.setZero();
                return;
            }
        }
    }
}

void LogdetFunctionWithConstraints::hessian(const Eigen::VectorXd &x, Eigen::MatrixXd &H)
{
    LogdetFunction::hessian(x, H);
    std::list<Eigen::MatrixXd> X = decondense(x);

    if(_chol.isPositive() && (_chol.vectorD().array() > 0).all()) {
        int s = 0, t = 0;
        for(const Eigen::MatrixXd &invXblock: _invXblocks) {
            int n = invXblock.rows();
            Eigen::VectorXd hij(n * n);
            for(int j = 0; j < n; j++) {
                for(int i = 0; i < n; i++, s++) {
                    singleVariableHessianBlock(invXblock, i, j, hij);
                    H.block(s, t, 1, n * n) += _rho * hij.transpose();
                }
            }
            t += n * n;
        }
    }
}

void LogdetFunctionWithConstraints::singleVariableHessianBlock(const Eigen::MatrixXd &P, int i, int j, Eigen::VectorXd &hij) const
{
    int s = 0;
    for(int v = 0; v < P.cols(); v++) {
        for(int u = 0; u < P.rows(); u++) {
            hij(s++) = P(u,i) * P(j,v);
        }
    }
}
