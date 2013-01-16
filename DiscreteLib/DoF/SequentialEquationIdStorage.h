/**
 * \file   SequentialEquationIdStorage.h
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief
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
#include <iostream>

#include "BaseLib/CodingTools.h"
#include "IEquationIdStorage.h"


namespace DiscreteLib
{

/**
 * \brief Storage assuming sequential and continuous index
 */
class SequentialEquationIdStorage : public IEquationIdStorage
{
public:

    /**
     * constructor
     * @param pt_id_start   Start point id
     * @param n             Size of total points
     */
    SequentialEquationIdStorage(std::size_t pt_id_start, std::size_t n)
    : _pt_id_start(pt_id_start), _n(n), _dof_begin(0), _delta_per_pt(1)
    {};

    ///
    virtual ~SequentialEquationIdStorage() {};

    ///
    virtual bool isSequential() const {return true;};

    ///
    virtual bool hasKey(std::size_t pt_id) const
    {
        return (_pt_id_start<=pt_id && pt_id<_pt_id_start+_n);
    }

    virtual bool hasValue(std::size_t eqs_id) const;

    virtual void key_range(std::size_t &i_start, std::size_t &i_end) const;

    virtual void activate(std::size_t pt_id, bool b);
    
    virtual bool isActive(std::size_t pt_id) const { return _deactive.count(pt_id)==0;};

    virtual void set(std::size_t /*pt_id*/, long /*eqs_id*/)
    {
        //invalid
        std::cout << "***Error: called invalid functions. SequentiallyMappedAddress::set()." << std::endl;
    }

    virtual std::size_t setAll(std::size_t address_start, std::size_t dn_pt);

    virtual std::size_t size() const {return _n;};

    virtual std::size_t address(std::size_t pt_id) const;

    virtual std::size_t key(std::size_t /*address_id*/) const
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
