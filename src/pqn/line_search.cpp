/*
 * line_search.cpp
 *
 *  Created on: Feb 14, 2014
 *      Author: mazuran
 */

#include "line_search.h"
#include "pqn_general.h"
#include <iostream>

int LineSearchSimpleBacktracking::findStep(
        Function &fun, PQNOptimizerOptions &options, const Eigen::VectorXd &d,
        double s,      double f, const Eigen::VectorXd &g, const Eigen::VectorXd &x,
        double &s_new, double &f_new,  Eigen::VectorXd &g_new,   Eigen::VectorXd &x_new)
{
    int funevals = 0;
    s_new = s;
    while(1) {
        x_new = x + s_new * d;
        double fold = fun.value(x);
        f_new = fun.value(x_new);

        if(s_new < 1e-12) {
            return -1;
        }

        funevals++;
        if(invalid(f_new) || f_new > f) {
            s_new /= 2;
        } else {
            fun.gradient(x_new, g_new);
            return funevals;
        }
    }
}

int LineSearchBacktracking::findStep(
        Function &fun, PQNOptimizerOptions &options, const Eigen::VectorXd &d,
        double s,      double f, const Eigen::VectorXd &g, const Eigen::VectorXd &x,
        double &s_new, double &f_new,  Eigen::VectorXd &g_new,   Eigen::VectorXd &x_new)
{
    /* TODO: Only geometric for now */
    int funevals = 0;
    const double decrease = 0.5, increase = 2.1;
    double t = decrease, g_directional = g.dot(d);
    s_new = s;
    while(1) {
        x_new = x + s_new * d;
        f_new = fun.value(x_new);
        funevals++;
        if(invalid(f_new) || f_new > f + options.sufficientDecrease * s_new * g_directional) {
            t = decrease;
        } else {
            fun.gradient(x_new, g_new);
            if(_cond == ARMIJO) return funevals;
            double g_directional_new = g_new.dot(d);
            if(g_directional_new < options.curvatureCondition * g_directional) {
                t = increase;
            } else {
                if(_cond == WOLFE) return funevals;
                /* if(cond == STRONG_WOLFE) */
                if(g_directional_new > -options.curvatureCondition * g_directional) {
                    t = decrease;
                } else {
                    return funevals;
                }
            }
        }
        s_new *= t;
    }
}

int LineSearchCubicBacktracking::findStep(
        Function &fun, PQNOptimizerOptions &options, const Eigen::VectorXd &d,
        double s,      double f, const Eigen::VectorXd &g, const Eigen::VectorXd &x,
        double &s_new, double &f_new,  Eigen::VectorXd &g_new,   Eigen::VectorXd &x_new)
{
    int funevals = 0;
    double gdotd = g.dot(d);
    x_new = x + s * d;
    f_new = fun.value(x_new);
    fun.gradient(x_new, g_new);
    funevals++;

    while(f_new > f + options.sufficientDecrease * s * gdotd || invalid(f_new)) {
        double backup = s;
        if(invalid(f_new) || invalid(g_new)) {
            s /= 2;
        } else {
            s = twoPointCubicInterpolation(
                    Eigen::Vector3d(0, f, gdotd), Eigen::Vector3d(s, f_new, g_new.dot(d)));
        }

        if(s < backup * 1e-3) {
            s = backup * 1e-3;
        } else if(s > backup * 0.6) {
            s = backup * 0.6;
        }

        if(s * d.cwiseAbs().sum() < options.optimalityTolerance) {
            if(options.verbose) std::cerr << "Line search failed" << std::endl;
            s = 0;
            f_new = f;
            g_new = g;
            break;
        }

        x_new = x + s * d;
        f_new = fun.value(x_new);
        fun.gradient(x_new, g_new);
        funevals++;
    }

    s_new = s;

    return funevals;
}

int LineSearchDirectionalHessian::findStep(
        Function &fun, PQNOptimizerOptions &options, const Eigen::VectorXd &d,
        double s,      double f, const Eigen::VectorXd &g, const Eigen::VectorXd &x,
        double &s_new, double &f_new,  Eigen::VectorXd &g_new,   Eigen::VectorXd &x_new)
{
    double h = fun.directionalHessian(x, d);
    s_new = - g.dot(d) / h;
    return _fallback.findStep(fun, options, d, s_new, f, g, x, s_new, f_new, g_new, x_new);
    /*
    x_new = x + s_new * d;
    f_new = fun.value(x_new);
    fun.gradient(x_new, g_new);

    if(invalid(f_new) || invalid(g_new)) {
        return 1 + _fallback.findStep(fun, options, d, s, f, g, x, s_new, f_new, g_new, x_new);
    } else {
        return 1;
    }
    */
}
