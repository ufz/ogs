/**
 * \file   SequentialEquationIdStorage.h
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

#ifndef SEQUENTIALEQUATIONIDSTORAGE_H_
#define SEQUENTIALEQUATIONIDSTORAGE_H_

#include <set>
#include <algorithm>
#include <cassert>
#include <iostream>

#include "BaseLib/CodingTools.h"
#include "IEquationIdStorage.h"


namespace DiscreteLib
{

/**
 *
 */
class SequentialEquationIdStorage : public IEquationIdStorage
{
public:

    /**
     *
     * @param pt_id_start
     * @param n
     */
    SequentialEquationIdStorage(std::size_t pt_id_start, std::size_t n)
    : _pt_id_start(pt_id_start), _n(n), _dof_begin(0), _delta_per_pt(1)
    {};

    virtual ~SequentialEquationIdStorage() {};

    bool isSequential() const {return true;};

    bool hasKey(std::size_t pt_id) const
    {
        return (_pt_id_start<=pt_id && pt_id<_pt_id_start+_n);
    }

    bool hasValue(std::size_t eqs_id) const;

    void key_range(std::size_t &i_start, std::size_t &i_end) const;

    void activate(std::size_t pt_id, bool b);
    
    bool isActive(std::size_t pt_id) const { return _deactive.count(pt_id)==0;};

    void set(std::size_t /*pt_id*/, long /*eqs_id*/)
    {
        //invalid
        std::cout << "***Error: called invalid functions. SequentiallyMappedAddress::set()." << std::endl;
    }

    std::size_t setAll(std::size_t address_start, std::size_t dn_pt);

    std::size_t size() const {return _n;};

    std::size_t address(std::size_t pt_id) const;

    std::size_t key(std::size_t /*address_id*/) const
    {
        //TODO
        return index_npos;
    }

private:
    std::size_t _pt_id_start;
    std::size_t _n;
    std::size_t _dof_begin;
    std::size_t _delta_per_pt;
    std::set<std::size_t> _deactive;
};

} //end

#endif //SEQUENTIALEQUATIONIDSTORAGE_H_
