/*
 * line_search.h
 *
 *  Created on: Feb 14, 2014
 *      Author: mazuran
 */

#ifndef LINE_SEARCH_H_
#define LINE_SEARCH_H_

#include <Eigen/Core>
#include "pqn_optimizer_options.h"
#include "function.h"

class LineSearch {
public:
    virtual ~LineSearch() {}
    virtual int findStep(

            Function &fun,
            PQNOptimizerOptions &options,
            const Eigen::VectorXd &d,

            double s,
            double f,
            const Eigen::VectorXd &g,
            const Eigen::VectorXd &x,

            double &s_new,
            double &f_new,
            Eigen::VectorXd &g_new,
            Eigen::VectorXd &x_new

    ) = 0;
};

class LineSearchSimpleBacktracking : public LineSearch {
public:
    int findStep(
            Function &fun, PQNOptimizerOptions &options, const Eigen::VectorXd &d,
            double s,      double f, const Eigen::VectorXd &g, const Eigen::VectorXd &x,
            double &s_new, double &f_new,  Eigen::VectorXd &g_new,   Eigen::VectorXd &x_new);
};


class LineSearchBacktracking : public LineSearch {
public:
    enum Condition     { ARMIJO, WOLFE, STRONG_WOLFE };
    enum Interpolation { GEOMETRIC, QUADRATIC, CUBIC };

    LineSearchBacktracking(Condition cond = WOLFE, Interpolation interp = CUBIC) :
        _cond(cond), _interp(interp)
    {}

    int findStep(
            Function &fun, PQNOptimizerOptions &options, const Eigen::VectorXd &d,
            double s,      double f, const Eigen::VectorXd &g, const Eigen::VectorXd &x,
            double &s_new, double &f_new,  Eigen::VectorXd &g_new,   Eigen::VectorXd &x_new);
private:
    Condition _cond;
    Interpolation _interp;
};

class LineSearchCubicBacktracking : public LineSearch {
public:
    int findStep(
            Function &fun, PQNOptimizerOptions &options, const Eigen::VectorXd &d,
            double s,      double f, const Eigen::VectorXd &g, const Eigen::VectorXd &x,
            double &s_new, double &f_new,  Eigen::VectorXd &g_new,   Eigen::VectorXd &x_new);
};

class LineSearchDirectionalHessian : public LineSearch {
public:
    LineSearchDirectionalHessian() : _fallback(LineSearchBacktracking::ARMIJO) {}

    int findStep(
            Function &fun, PQNOptimizerOptions &options, const Eigen::VectorXd &d,
            double s,      double f, const Eigen::VectorXd &g, const Eigen::VectorXd &x,
            double &s_new, double &f_new,  Eigen::VectorXd &g_new,   Eigen::VectorXd &x_new);
private:
    LineSearchBacktracking _fallback;
};

#endif /* LINE_SEARCH_H_ */
