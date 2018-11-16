/*
 * pqn_optimizer.cpp
 *
 *  Created on: 24/gen/2014
 *      Author: Mladen Mazuran
 */

#include "pqn_optimizer.h"
#include "pqn_general.h"
#include <stdio.h>
#include <iostream>
#include <Eigen/Cholesky>

PQNOptimizer::PQNOptimizer() : _search(new LineSearchBacktracking)
{
}

PQNOptimizer::~PQNOptimizer()
{
    delete _search;
}

void PQNOptimizer::setLineSearch(LineSearch *search)
{
    delete _search;
    _search = search;
}

double PQNOptimizer::optimize(Function &fun, Eigen::VectorXd &x)
{
    Eigen::VectorXd g(x.rows()), g_old, x_old, precon(x.rows());
    double f = fun.value(x), f_old = 0;
    fun.gradient(x, g);
    int i = 0;
    _funevals = 1;

    if(_opts.verbose) {
        fprintf(stderr, "%10s %10s %15s %15s %15s\n", "Iteration", "FunEvals", "Step Length", "Function Val", "Opt Cond");
        fprintf(stderr, "%10d %10d %15s %15.8f %15s\n", i++, _funevals, " ", f, " ");
    }

    _S.clear();
    _Y.clear();

    while(i < _opts.maxIters || _opts.maxIters == 0) {
        double step;
        Eigen::VectorXd d;
        if(_opts.useHessian) {
            Eigen::MatrixXd H(g.rows(), g.rows());
            fun.hessian(x, H);
            // TODO: Assumes PSD, which is true for convex functions, but, yeah
            Eigen::LLT<Eigen::MatrixXd> chol(H);
            d = chol.solve(-g);
        } else if(i == 1) {
            d = -g;
            _hdiag = 1;
        } else {
            lbfgsUpdate(g - g_old, x - x_old);
            precon.setConstant(precon.rows(), _hdiag);
            if(_opts.usePreconditioner) {
                fun.preconditioner(x, precon);
            }
            d = lbfgsSearchDirection(-g, precon);
        }

        double gdotd = g.dot(d);
        if(std::abs(gdotd) < _opts.optimalityTolerance) {
            if(_opts.verbose) std::cerr << "Directional derivative below optimality tolerance" << std::endl;
            return f;
        }

        /*
        if(i == 0) {
            step = 1 / g.cwiseAbs().sum();
        } else {
            step = std::min(1., 2 * (f - f_old) / gdotd);
        }
        */
        step = 1;

        f_old = f;
        g_old = g;
        x_old = x;

        double f_new = f;
        Eigen::VectorXd x_new(x.rows()), g_new(g.rows());
        int ret = _search->findStep(fun, _opts, d, step, f, g, x, step, f_new, g_new, x_new);

        if(ret < 0) {
            if(_opts.verbose) {
                std::cerr << "Line search failed" << std::endl;
            }
            return INFINITY;
        }

        _funevals += ret;

        double optcond = g.cwiseAbs().sum(); // TODO: Check

        if(_opts.verbose) {
            fprintf(stderr, "%10d %10d %15.8f %15.8f %15.8f\n", i, _funevals, step, f, optcond);
        }

        x = x_new;
        f = f_new;
        g = g_new;

        if(optcond < _opts.optimalityTolerance) {
            if(_opts.verbose) std::cerr << "First-order optimality condition below optimality tolerance" << std::endl;
            return f;
        }

        if(step * d.cwiseAbs().sum() < _opts.optimalityTolerance) {
            if(_opts.verbose) std::cerr << "Step size below optimality tolerance" << std::endl;
            return f;
        }

        //if(std::abs(f - f_old) / std::abs(f_old) < _opts.optimalityTolerance) {
        if(std::abs(f - f_old) < _opts.optimalityTolerance) {
            if(_opts.verbose) std::cerr << "Function value changing by less than optimality tolerance" << std::endl;
            return f;
        }

        i++;
    }

    return f;
}

void PQNOptimizer::lbfgsUpdate(const Eigen::VectorXd &y, const Eigen::VectorXd &s)
{
    if(y.dot(s) > 1e-10) {
        if(_Y.size() >= (size_t) _opts.lbfgsCorrections) {
            _Y.pop_front();
            _S.pop_front();
        }
        _Y.push_back(y);
        _S.push_back(s);
        _hdiag = y.dot(s) / y.dot(y);
    } /* TODO: ??? */
}

Eigen::VectorXd PQNOptimizer::lbfgsSearchDirection(const Eigen::VectorXd &g, const Eigen::VectorXd &precon)
{
    Eigen::VectorXd rho(_Y.size());
    for(size_t i = 0; i < _Y.size(); i++) {
        rho[i] = 1 / _Y[i].dot(_S[i]);
    }

    Eigen::VectorXd q = g, alpha(_Y.size());

    for(int i = int(_Y.size()) - 1; i >= 0; i--) {
        alpha[i] = rho[i] * _S[i].dot(q);
        q -= alpha[i] * _Y[i];
    }

    Eigen::VectorXd r = precon.cwiseProduct(q);

    for(size_t i = 0; i < _Y.size(); i++) {
        double beta = rho[i] * _Y[i].dot(r);
        r += (alpha[i] - beta) * _S[i];
    }

    return r;
}
