/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file OrderedMap.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#ifndef ORDERED_MAP_H_
#define ORDERED_MAP_H_

#include <map>
#include <vector>

namespace BaseLib
{

/**
 * \brief Map container keeping the inserted order of values
 *
 * \tparam T_KEY    Type of key
 * \tparam T_VAL    Type of value
 */
template <typename T_KEY, typename T_VAL>
class OrderedMap
{
public:
    typedef typename std::map<T_KEY, T_VAL>::iterator iterator;
    typedef typename std::map<T_KEY, T_VAL>::const_iterator const_iterator;

    ///
    OrderedMap() {};

    /// 
    ~OrderedMap() {};

    /**
     * insert a key and value
     * 
     * \param key   the key
     * \param val   corresponding value
     */
    void insert(const T_KEY &key, const T_VAL &val)
    {
        if (_map.count(key)==0)
            _vec.push_back(key);
        _map[key] = val;
    }

    /**
     * return the number of entries matched with the given key
     *
     * \param val   corresponding value
     * \return  the numbe of entries matched
     */
    size_t count(const T_KEY &key) const
    {
        return _map.count(key);
    }

    /**
     * get the number of entries
     *
     * \return  the numbe of entries
     */
    size_t size() const {return _vec.size();};

    /// return the beginning of this container
    iterator begin() {return _map.begin();};
    /// return the end of this container
    iterator end() {return _map.end();};

    /// return the beginning of this container
    const_iterator begin() const {return _map.begin();};
    /// return the end of this container
    const_iterator end() const {return _map.end();};

    /// clear entries
    void clear()
    {
        _map.clear();
        _vec.clear();
    }

    /// find an iterator to an entry matched with the given key
	iterator find(const T_KEY & key) 
	{
        return _map.find(key);
    }

    /// find an iterator to an entry matched with the given key
	const_iterator find(const T_KEY & key) const
	{
        return _map.find(key);
    }

    /// access an iterator by index
	iterator operator[] (size_t idx)
    {
        return _map.find(_vec.at(idx));
    }

    /// access an iterator by index
	const_iterator operator[] (size_t idx) const
    {
        return _map.find(_vec.at(idx));
    }

private:
    std::map<T_KEY, T_VAL> _map;
    std::vector<T_KEY> _vec;
};

}

#endif // end ORDERED_MAP_H_