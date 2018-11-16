/*
 * topology_provider_base.h
 *
 *  Created on: Oct 16, 2014
 *      Author: mazuran
 */

#ifndef TOPOLOGY_PROVIDER_BASE_H_
#define TOPOLOGY_PROVIDER_BASE_H_

#include "topology_provider.h"
#include "utils.h"

class TopologyProviderBase : public TopologyProvider
{
public:
    TopologyProviderBase() : TopologyProvider() {}
    virtual ~TopologyProviderBase() {}

    virtual bool applicable(
            const VertexInfoSet &vinfo,
            const EdgeInfoSet &einfo) {
        return setDifference(vinfo, _okvertices).size() == 0 &&
                setDifference(einfo, _okedges).size() == 0;
    }

    virtual bool requiresOptimization() const {
        return true;
    }

protected:
    VertexInfoSet _okvertices;
    EdgeInfoSet _okedges;
};

#endif /* TOPOLOGY_PROVIDER_BASE_H_ */
