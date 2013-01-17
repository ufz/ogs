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
 * \brief Equation ID storage with sequential point and equation IDs
 *
 * This storage assumes that equation ID can be analytically determined for
 * a given equation ID. Therefore it can save memory and work faster.
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

    /**
     *
     * @param pt_id
     * @return
     */
    virtual bool hasPoint(std::size_t pt_id) const
    {
        return (_pt_id_start<=pt_id && pt_id<_pt_id_start+_n);
    }

    /**
     *
     * @param eqs_id
     * @return
     */
    virtual bool hasEquationID(std::size_t eqs_id) const;

    /**
     *
     * @param i_start
     * @param i_end
     */
    virtual void getPointRange(std::size_t &i_start, std::size_t &i_end) const;

    /**
     *
     * @param pt_id
     * @param b
     */
    virtual void activate(std::size_t pt_id, bool b);
    
    /**
     *
     * @param pt_id
     * @return
     */
    virtual bool isActive(std::size_t pt_id) const { return _deactive.count(pt_id)==0;};

    /**
     *
     * @param
     * @param
     */
    virtual void set(std::size_t /*pt_id*/, std::size_t /*eqs_id*/)
    {
        //nothing happen
    }

    /**
     * renumbering all equation index
     * @param eqs_id_start
     * @param delta_pt
     * @return
     */
    virtual std::size_t setAll(std::size_t eqs_id_start, std::size_t delta_pt);

    /**
     *
     * @return
     */
    virtual std::size_t size() const {return _n;};

    /**
     *
     * @param pt_id
     * @return
     */
    virtual std::size_t equationID(std::size_t pt_id) const;

    /**
     *
     * @param
     * @return
     */
    virtual std::size_t pointID(std::size_t /*address_id*/) const
    {
        //TODO
        return index_npos;
    }

    /// print debug info
    virtual void printout() const
    {
    }

private:
    const std::size_t _pt_id_start;
    const std::size_t _n;
    std::size_t _dof_begin;
    std::size_t _delta_per_pt;
    std::set<std::size_t> _deactive;
};

} //end

#endif //SEQUENTIALEQUATIONIDSTORAGE_H_
