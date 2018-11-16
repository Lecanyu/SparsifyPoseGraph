/*
 * logdet_function.h
 *
 *  Created on: May 22, 2014
 *      Author: mazuran
 */

#ifndef LOGDET_FUNCTION_H_
#define LOGDET_FUNCTION_H_

#include "pqn/pqn_optimizer.h"
#include <set>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include "optimizer.h"

class LogdetFunction : public Function
{
public:
    LogdetFunction(
            const JacobianMapping &mapping,
            const Eigen::MatrixXd &targetInformation);

    bool hasClosedFormSolution() const;
    std::set<int> chooseDimensions(const Eigen::MatrixXd &candidates) const;

    virtual double value(const Eigen::VectorXd &x);
    virtual void gradient(const Eigen::VectorXd &x, Eigen::VectorXd &g);

    virtual void hessian(const Eigen::VectorXd &x, Eigen::MatrixXd &H);

    Eigen::VectorXd educatedGuess();
    std::list<Eigen::MatrixXd> closedFormSolution() const;
    std::list<Eigen::MatrixXd> decondense(const Eigen::VectorXd &x) const;
    Eigen::VectorXd condense(const std::list<Eigen::MatrixXd> &X) const;

    Eigen::MatrixXd informationProduct(const std::list<Eigen::MatrixXd> &X) const;
    Eigen::MatrixXd informationProduct(const Eigen::VectorXd &x) const;
    Eigen::SparseMatrix<double> sparseJacobian() const;

protected:
    void singleVariableHessian(const Eigen::MatrixXd &P, int i, int j, Eigen::VectorXd &hij) const;

protected:
    const JacobianMapping &_mapping;
    const Eigen::MatrixXd &_target;
    double _logdet, _jacsize;
    Eigen::VectorXd _S;
    Eigen::MatrixXd _U, _xinv;
    Eigen::LDLT<Eigen::MatrixXd> _chol;
};

class LogdetFunctionWithConstraints : public LogdetFunction
{
public:
    LogdetFunctionWithConstraints(
            const JacobianMapping &mapping,
            const Eigen::MatrixXd &targetInformation);

    double rho() const { return _rho; }
    void setRho(double rho) { _rho = rho; }

    virtual double value(const Eigen::VectorXd &x);
    virtual void gradient(const Eigen::VectorXd &x, Eigen::VectorXd &g);
    virtual void hessian(const Eigen::VectorXd &x, Eigen::MatrixXd &H);

private:
    void singleVariableHessianBlock(const Eigen::MatrixXd &P, int i, int j, Eigen::VectorXd &hij) const;

private:
    double _rho;
    std::list<Eigen::MatrixXd> _invXblocks;
};

#endif /* LOGDET_FUNCTION_H_ */
