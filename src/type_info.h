/*
 * type_info.h
 *
 *  Created on: Oct 14, 2014
 *      Author: mazuran
 */

#ifndef TYPE_INFO_H_
#define TYPE_INFO_H_

#include <g2o/core/optimizable_graph.h>
#include <unordered_map>
#include <vector>
#include <Eigen/Core>
#include "vertex_info.h"
#include "edge_info.h"
#include <typeinfo>

class TypeInfo {
public:
    TypeInfo();
    virtual ~TypeInfo();
    
    static TypeInfo *instance();
    
    static void registerType(VertexInfo *vi);
    static void registerType(EdgeInfo *ei);

    const VertexInfo *operator()(const g2o::HyperGraph::Vertex *v) const;
    const VertexInfo *vertex(const std::type_info &t) const;
    const VertexInfo *vertex(const g2o::HyperGraph::Vertex *v) const;
    const VertexInfo *vertex(size_t id) const;

    const EdgeInfo *operator()(const g2o::HyperGraph::Edge *e) const;
    const EdgeInfo *edge(const std::type_info &t) const;
    const EdgeInfo *edge(const g2o::HyperGraph::Edge *e) const;
    const EdgeInfo *edge(size_t id) const;

    const UnknownVertexInfo *vertexUnknown() const { return _vunk; }
    const UnknownEdgeInfo   *edgeUnknown()   const { return _eunk; }

private:
    std::unordered_map<size_t, VertexInfo *> _vinfos;
    std::unordered_map<size_t, EdgeInfo *>   _einfos;
    UnknownVertexInfo        *_vunk;
    UnknownEdgeInfo          *_eunk;
    static TypeInfo *_instance;
};

#define typeinfo (*TypeInfo::instance())

template <typename T>
struct TypeInfoRegisterProxy {
    TypeInfoRegisterProxy() {
        TypeInfo::registerType(new T);
    }
};

struct TypeInfoLinkProxy {
    TypeInfoLinkProxy(void (*f)()) {
        f();
    }
};

#define G2S_REGISTER_TYPE(cls)                                  \
namespace __internal {                                          \
    void g2s_link_##cls() {}                                    \
    TypeInfoRegisterProxy<cls> g2s_proxy_##cls;                 \
}

#define G2S_USE_TYPE(cls)                                       \
namespace __internal {                                          \
    extern void g2s_link_##cls();                               \
    static TypeInfoLinkProxy g2s_forcer_##cls(g2s_link_##cls);  \
}

#endif /* TYPE_INFO_H_ */
