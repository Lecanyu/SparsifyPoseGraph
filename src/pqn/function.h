/*
 * function.h
 *
 *  Created on: Feb 14, 2014
 *      Author: mazuran
 */

#ifndef FUNCTION_H_
#define FUNCTION_H_

class Function {
public:
    /* Returns the value of the function evaluated at point x */
    virtual double value(const Eigen::VectorXd &x) = 0;
    /* Fills g with the gradient of the function evaluated at point x */
    virtual void gradient(const Eigen::VectorXd &x, Eigen::VectorXd &g) = 0;

    virtual void hessian(const Eigen::VectorXd &x, Eigen::MatrixXd &H) {}

    /* Fills prec with a diagonal preconditioner for the function evaluated at x (optional) */
    virtual void preconditioner(const Eigen::VectorXd &x, Eigen::VectorXd &prec) {}
    /* Fills projected with the feasible argument closest to x in terms of Euclidean norm (optional) */
    virtual bool project(const Eigen::VectorXd &x, Eigen::VectorXd &projected) { projected = x; return true; }
    /* Returns the Hessian along the direction d evaluated at point x (optional) */
    virtual double directionalHessian(const Eigen::VectorXd &x, const Eigen::VectorXd &d) { return 1; }
};

#endif /* FUNCTION_H_ */
