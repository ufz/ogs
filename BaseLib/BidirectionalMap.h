/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file BidirectionalMap.h
 *
 * Created on 2012-02-24 by Norihiro Watanabe
 */

#ifndef BIDIRECTIONAL_MAP_H_
#define BIDIRECTIONAL_MAP_H_

#include <map>
#include <stdexcept>

namespace BaseLib
{

/**
 * \brief A bidrectional map class for two sets of unique keys A and B
 *
 * \tparam T1   Value type of the first key
 * \tparam T2   Value type of the second key
 */
template <typename T1, typename T2>
class BidirectionalMap
{
public:
    typedef std::map<T1, T2> MapA;
    typedef std::map<T2, T1> MapB;

    /**
     * Create a new container
     */
    BidirectionalMap() {};

    /**
     * Create a new container and copy the existing container
     *
     * \param src   source object
     */
    explicit BidirectionalMap(const BidirectionalMap &src)
    : _map1(src._map1), _map2(src._map2)
    {};

    /**
     * Copy the given container to this
     *
     * \param src   source object
     */
    BidirectionalMap &operator=(const BidirectionalMap &src)
    {
        _map1 = src._map1;
        _map2 = src._map2;
        return *this;
    }

    /**
     * insert a set of keys
     *
     * \param v1    Key1
     * \param v2    Key2
     */
    void insert(const T1 &v1, const T2 &v2)
    {
        _map1[v1] = v2;
        _map2[v2] = v1;
    }

    /**
     * return the number of elements found with key1
     *
     * \param v    Key1
     */
    std::size_t countInA(const T1 &v) const { return _map1.count(v);};

    /**
     * return the number of elements found with key2
     *
     * \param v    Key2
     */
    std::size_t countInB(const T2 &v) const { return _map2.count(v);};

    /**
     * return size of entries
     */
    std::size_t size() const { return _map1.size(); };

    /**
     * clear this container
     */
    void clear() 
    {
        _map1.clear(); 
        _map2.clear(); 
    };

    /**
     * return key2 corresponding to the given key1
     *
     * This function throws an exception if the key1 is not found.
     * \param v    Key1
     */
    const T2 mapAtoB(const T1 &v) const 
    {
        typename MapA::const_iterator itr = _map1.find(v);
        if (itr!=_map1.end())
            return itr->second;
        else
            throw std::out_of_range("the given key for A is not found.");
    };

    /**
     * return key1 corresponding to the given key2
     *
     * This function throws an exception if the key1 is not found.
     * \param v    Key2
     */
    const T1 mapBtoA(const T2 &v) const 
    {
        typename MapB::const_iterator itr = _map2.find(v);
        if (itr!=_map2.end())
            return itr->second;
        else
            throw std::out_of_range("the given key for B is not found.");
    };

private:
    MapA _map1;
    MapB _map2;
};


}

#endif // end BIDIRECTIONAL_MAP_H_
