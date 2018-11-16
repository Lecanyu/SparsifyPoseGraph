/*
 * type_info.cpp
 *
 *  Created on: Oct 14, 2014
 *      Author: mazuran
 */

#include "type_info.h"

TypeInfo *TypeInfo::_instance = NULL;

TypeInfo::TypeInfo() : _vunk(new UnknownVertexInfo), _eunk(new UnknownEdgeInfo)
{
#ifdef G2S_VERBOSE
    std::cerr << "Successfully called TypeInfo constructor" << std::endl;
#endif
}

TypeInfo::~TypeInfo()
{
    for(std::unordered_map<size_t, VertexInfo *>::iterator it = _vinfos.begin();
            it != _vinfos.end(); ++it) {
        delete it->second;
    }
    for(std::unordered_map<size_t, EdgeInfo *>::iterator it = _einfos.begin();
            it != _einfos.end(); ++it) {
        delete it->second;
    }
    delete _vunk;
    delete _eunk;
}

TypeInfo *TypeInfo::instance()
{
    if(!_instance) {
        _instance = new TypeInfo;
    }
    return _instance;
}

void TypeInfo::registerType(VertexInfo *info)
{
    instance()->_vinfos.insert(std::make_pair(info->id(), info));
#ifdef G2S_VERBOSE
    std::cerr << "Registered vertex type " << info->name() << " with id " << info->id() << std::endl;
#endif
}

void TypeInfo::registerType(EdgeInfo *info)
{
    instance()->_einfos.insert(std::make_pair(info->id(), info));
#ifdef G2S_VERBOSE
    std::cerr << "Registered edge type " << info->name() << " with id " << info->id() << std::endl;
#endif
}

const VertexInfo *TypeInfo::operator()(const g2o::HyperGraph::Vertex *v) const
{
    return vertex(v);
}

const VertexInfo *TypeInfo::vertex(const g2o::HyperGraph::Vertex *v) const
{
    return vertex(typeid(*v).hash_code());
}

const VertexInfo *TypeInfo::vertex(size_t id) const
{
    std::unordered_map<size_t, VertexInfo *>::const_iterator it = _vinfos.find(id);
    if (it != _vinfos.end()) return it->second;
    else return _vunk;
}

const VertexInfo *TypeInfo::vertex(const std::type_info &t) const
{
    return vertex(t.hash_code());
}

const EdgeInfo *TypeInfo::operator()(const g2o::HyperGraph::Edge *e) const
{
    return edge(e);
}

const EdgeInfo *TypeInfo::edge(const g2o::HyperGraph::Edge *e) const
{
    return edge(typeid(*e).hash_code());
}

const EdgeInfo *TypeInfo::edge(size_t id) const
{
    std::unordered_map<size_t, EdgeInfo *>::const_iterator it = _einfos.find(id);
    if (it != _einfos.end()) return it->second;
    else return _eunk;
}

const EdgeInfo *TypeInfo::edge(const std::type_info &t) const
{
    return edge(t.hash_code());
}

