/*
 * decimation.h
 *
 *  Created on: Oct 23, 2014
 *      Author: mazuran
 */

#ifndef DECIMATION_H_
#define DECIMATION_H_

#include <vector>

struct DecimateOptions {
    int sparsity;
    int clusterSize;
};

std::vector<int> clusterDecimate(int last, int endvert, const DecimateOptions &opts);
std::vector<int> onlineDecimate(int last, int endvert, const DecimateOptions &opts);
std::vector<int> globalDecimate(int lastid, int endvert, const DecimateOptions &opts);

typedef std::vector<int> (*DecimateFunction)(int last, int endvert, const DecimateOptions &opts);

#endif /* DECIMATION_H_ */
