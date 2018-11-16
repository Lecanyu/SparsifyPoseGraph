/*
 * pqn_optimizer.h
 *
 *  Created on: 24/gen/2014
 *      Author: Mladen Mazuran
 */

#ifndef PQN_OPTIMIZER_H_
#define PQN_OPTIMIZER_H_

#include <Eigen/Core>
#include <deque>
#include "pqn_optimizer_options.h"
#include "line_search.h"
#include "function.h"

class PQNOptimizer
{
public:
    PQNOptimizer();
    virtual ~PQNOptimizer();

    double optimize(Function &fun, Eigen::VectorXd &x);
    void lbfgsUpdate(const Eigen::VectorXd &y, const Eigen::VectorXd &s);

    Eigen::VectorXd lbfgsSearchDirection(const Eigen::VectorXd &g, const Eigen::VectorXd &precon);

    void setOptions(const PQNOptimizerOptions &opts) { _opts = opts; }
    const PQNOptimizerOptions &options() const { return _opts; }

    void setLineSearch(LineSearch *search);
    LineSearch *lineSearch() { return _search; }
    const LineSearch *lineSearch() const { return _search; }

private:
    LineSearch *_search;
    PQNOptimizerOptions _opts;
    int _funevals;
    double _hdiag;
    std::deque<Eigen::VectorXd> _Y, _S;
};

#endif /* PQNOPTIMIZER_H_ */
