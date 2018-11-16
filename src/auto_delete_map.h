/*
 * auto_delete_map.h
 *
 *  Created on: 11/gen/2014
 *      Author: Mladen Mazuran
 */

#ifndef AUTO_DELETE_MAP_H_
#define AUTO_DELETE_MAP_H_

#include <map>

template <typename K, typename T>
class AutoDeleteMap : public std::map<K, T>
{
public:
    typedef typename std::map<K, T>::iterator iterator;
    typedef typename std::map<K, T>::size_type size_type;

    virtual ~AutoDeleteMap() {
        for(auto v: *this) {
            delete v.second;
        }
    }

    void erase(iterator which) {
        delete which->second;
        std::map<K, T>::erase(which);
    }

    size_type erase(const K &which) {
        delete (*this)[which];
        return std::map<K, T>::erase(which);
    }

    void erase(iterator begin, iterator end) {
        iterator it = begin;
        while(it != end) {
            delete it->second;
            ++it;
        }
        std::map<K, T>::erase(begin, end);
    }

};


#endif /* AUTO_DELETE_MAP_H_ */
