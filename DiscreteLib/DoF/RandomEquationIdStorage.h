/**
 * \file   RandomEquationIdStorage.h
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief  Helper macros.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#ifndef RANDOMEQUATIONIDSTORAGE_H_
#define RANDOMEQUATIONIDSTORAGE_H_

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cassert>
#include <iostream>

#include <boost/bimap/bimap.hpp>

#include "BaseLib/CodingTools.h"
#include "IEquationIdStorage.h"

namespace DiscreteLib
{

/**
 *
 */
class RandomEquationIdStorage : public IEquationIdStorage
{
public:
    /**
     *
     * @param sorted_list_pt_id
     */
    RandomEquationIdStorage(const std::vector<std::size_t> &sorted_list_pt_id);

    /**
     *
     * @param pt_id_start
     * @param n
     */
    RandomEquationIdStorage(std::size_t pt_id_start, std::size_t n);

    ///
    virtual ~RandomEquationIdStorage() {};

    bool isSequential() const {return false;};

    bool hasKey(std::size_t pt_id) const;

    bool hasValue(std::size_t eqs_id) const;

    void key_range(std::size_t &i_start, std::size_t &i_end) const;

    void activate(std::size_t pt_id, bool b);

    bool isActive(std::size_t pt_id) const { return _deactive.count(pt_id)==0;};

    void set(std::size_t pt_id, long eqs_id);

    std::size_t setAll(std::size_t dof_start, std::size_t delta_per_pt);

    std::size_t size() const {return _map_pt2eqs.size();};

    std::size_t address(std::size_t pt_id) const;

    std::size_t key(std::size_t address_id) const;

private:
    typedef boost::bimaps::bimap<std::size_t, std::size_t> bimap_t;
    typedef bimap_t::value_type value_type;
    bimap_t _map_pt2eqs;
    std::vector<std::size_t> _list_pt_id;
    std::set<std::size_t> _deactive;
};

} //end


#endif //RANDOMEQUATIONIDSTORAGE_H_
