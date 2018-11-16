/*
 * evaluate.h
 *
 *  Created on: Oct 23, 2014
 *      Author: mazuran
 */

#ifndef EVALUATE_H_
#define EVALUATE_H_

#include "decimation.h"
#include "graph_wrapper.h"
#include <string>
#include <vector>

struct EvaluateInfo
{
    enum Algorithm {
        NFR, GLC, None
    };

    std::string g2oname;
    std::string destdir;
    Algorithm algorithm;
    bool useChi2;
    int kldPeriod;
    DecimateFunction decimate;
    DecimateOptions decimateOptions;
    SparsityOptions sparsityOptions;
};

void evaluate(GraphWrapper *gw, const EvaluateInfo &info);

void parallelEvaluate(
        const std::vector<EvaluateInfo> &joblist,
        int threads = 0);

#endif /* EVALUATE_H_ */
