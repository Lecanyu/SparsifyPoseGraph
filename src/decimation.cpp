/*
 * decimation.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: mazuran
 */

#include "decimation.h"
#include <cmath>

std::vector<int> clusterDecimate(int last, int endvert, const DecimateOptions &opts)
{
    if(((last - 4) % opts.clusterSize == 0 && last > 4) || last == endvert) {
        std::vector<int> ret;
        for(int i = int(std::ceil((last - 5) / (double) opts.clusterSize) -
                1) * opts.clusterSize + 5; i <= last; i++) {
            if(i % opts.sparsity > 0) {
                ret.push_back(i);
            }
        }
        return ret;
    } else {
        return std::vector<int>();
    }
}

std::vector<int> onlineDecimate(int last, int, const DecimateOptions &opts)
{
    if(last % opts.sparsity == 0) {
        return std::vector<int>();
    } else {
        return std::vector<int>({ last });
    }
}

std::vector<int> globalDecimate(int last, int endvert, const DecimateOptions &opts)
{
    if(last == endvert) {
        std::vector<int> which;
        for(int i = 4; i <= endvert; i++) {
            if(i % opts.sparsity != 0) {
                which.push_back(i);
            }
        }
        return which;
    } else {
        return std::vector<int>();
    }
}
