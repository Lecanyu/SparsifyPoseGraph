/*
 * kld_compare.cpp
 *
 *  Created on: Jan 20, 2014
 *      Author: mazuran
 */

#include <iostream>
#include "graph_wrapper_g2o.h"

int main(int argc, char **argv)
{
    if(argc < 3) {
        return -1;
    } else {
        GraphWrapperG2O *baseline       = new GraphWrapperG2O(argv[1], false);
        GraphWrapperG2O *incremental    = new GraphWrapperG2O(argv[2], false);
        
        incremental->optimize();
        baseline->optimize();

        incremental->push();
        baseline->push();

//        std::cout << baseline->estimateDifference(incremental) << std::endl << std::endl;
//        std::cout << baseline->estimateDifference(incremental).cwiseAbs().maxCoeff() << std::endl << std::endl;

        std::cout << "adapt marginalized to incremental" << std::endl;
        std::cout << "chi2 before = " << baseline->chi2() << std::endl;
        std::cout << "kld before  = " << baseline->kullbackLeibler(incremental) << std::endl;
        baseline->optimizeWithReference(incremental);
        std::cout << "chi2 after  = " << baseline->chi2() << std::endl;
        std::cout << "kld after   = " << baseline->kullbackLeibler(incremental) << std::endl << std::endl;
        
        incremental->pop();
        baseline->pop();

        std::cout << "adapt incremental to marginalized" << std::endl;
        std::cout << "chi2 before = " << incremental->chi2() << std::endl;
        std::cout << "kld before  = " << baseline->kullbackLeibler(incremental) << std::endl;
        incremental->optimizeWithReference(baseline);
        std::cout << "chi2 after  = " << incremental->chi2() << std::endl;
        std::cout << "kld after   = " << baseline->kullbackLeibler(incremental) << std::endl;
        
        delete baseline;
        delete incremental;
        return 0;
    }
}
