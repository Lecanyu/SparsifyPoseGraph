/*
 * sparsity_options.h
 *
 *  Created on: Oct 16, 2014
 *      Author: mazuran
 */

#ifndef SPARSITY_OPTIONS_H_
#define SPARSITY_OPTIONS_H_

struct SparsityOptions {
    enum SparsityTopology {
        Tree, Subgraph, CliqueySubgraph, Dense, CliqueyDense
    };

    enum LinearizationPoint {
        Local, Global
    };

    SparsityTopology topology;
    double chordRatio;
    LinearizationPoint linPoint;
    bool includeIntraClique;

    SparsityOptions() :
        topology(Tree),
        chordRatio(1),
        linPoint(Local),
        includeIntraClique(true) {}
};

#endif /* SPARSITY_OPTIONS_H_ */
