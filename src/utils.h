#ifndef MATRIXUTILS_H_
#define MATRIXUTILS_H_

#include <g2o/core/sparse_block_matrix.h>
#include <g2o/core/sparse_optimizer.h>
#include <isam/Slam.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <map>
#include <set>
#include <list>
#include <vector>

Eigen::VectorXd selectVariables(const Eigen::VectorXd &original, const std::vector<int> &variables);
Eigen::MatrixXd selectVariables(const Eigen::MatrixXd &original, const std::vector<int> &variables);
Eigen::MatrixXd selectVariables(const Eigen::MatrixXd &original, const std::vector<int> &variables1, const std::vector<int> &variables2);
Eigen::SparseMatrix<double> selectVariables(const Eigen::SparseMatrix<double> &original,
        const std::vector<int> &variables);
Eigen::SparseMatrix<double> selectVariables(const Eigen::SparseMatrix<double> &original,
        const std::vector<int> &variables1, const std::vector<int> &variables2);

//Eigen::MatrixXd computeDiff(const Eigen::VectorXd &mx, const Eigen::VectorXd &my, bool is2d = true, bool manifold = false);

enum KullbackMode {
    InformationInformation, InformationCovariance
};

double kullbackLeiblerDivergence(
        const Eigen::VectorXd &diff, const Eigen::MatrixXd &infox, const Eigen::MatrixXd &maty,
        KullbackMode mode = InformationCovariance);

std::vector<Eigen::Triplet<double> > triplets(const g2o::SparseBlockMatrix<Eigen::MatrixXd> &mat);
std::vector<Eigen::Triplet<double> > triplets(const Eigen::SparseMatrix<double> &mat);
Eigen::SparseMatrix<double> sparse(const isam::SparseSystem &mat);
Eigen::SparseMatrix<double> sparse(const g2o::SparseBlockMatrix<Eigen::MatrixXd> &mat);
Eigen::MatrixXd dense(const g2o::SparseBlockMatrix<Eigen::MatrixXd> &mat);

int fullSize(g2o::SparseOptimizer *so);
int fullSize(isam::Slam *isam);

Eigen::VectorXd stack(g2o::SparseOptimizer *so);
Eigen::VectorXd stack(isam::Slam *isam);

template <typename T, typename S>
std::set<T, S> setDifference(const std::set<T, S> &s1, const std::set<T, S> &s2)
{
    std::set<T, S> ret(s1.value_comp());
    for(typename std::set<T, S>::const_iterator it = s1.begin(); it != s1.end(); ++it) {
        if(s2.count(*it) == 0) {
            ret.insert(*it);
        }
    }
    return ret;
}

template <typename T, typename S>
std::set<T, S> setIntersection(const std::set<T, S> &s1, const std::set<T, S> &s2)
{
    std::set<T, S> ret(s1.value_comp());
    for(typename std::set<T, S>::const_iterator it = s1.begin(); it != s1.end(); ++it) {
        if(s2.count(*it) > 0) {
            ret.insert(*it);
        }
    }
    return ret;
}

template <typename Iter>
std::ostream &dumpIterable(std::ostream &s, Iter from, Iter to)
{
    s << "{";
    if(from != to) {
        s << *from;
        ++from;
    }
    while(from != to) {
        s << "," << *from;
        ++from;
    }
    return s << "}";
}

template <typename T, typename S>
std::ostream &operator<<(std::ostream &s, const std::pair<T,S> &p);
template <typename T>
std::ostream &operator<<(std::ostream &s, const std::list<T> &l);
template <typename T, typename S>
std::ostream &operator<<(std::ostream &s, const std::set<T, S> &l);
template <typename T>
std::ostream &operator<<(std::ostream &s, const std::vector<T> &l);

template <typename T, typename S>
std::ostream &operator<<(std::ostream &s, const std::pair<T,S> &p)
{
    return s << "(" << p.first << "," << p.second << ")";
}

template <typename T>
std::ostream &operator<<(std::ostream &s, const std::list<T> &l)
{
    return dumpIterable(s, l.begin(), l.end());
}

template <typename T, typename S>
std::ostream &operator<<(std::ostream &s, const std::set<T, S> &l)
{
    return dumpIterable(s, l.begin(), l.end());
}

template <typename T>
std::ostream &operator<<(std::ostream &s, const std::vector<T> &l)
{
    return dumpIterable(s, l.begin(), l.end());
}


#endif /* MATRIXUTILS_H_ */
