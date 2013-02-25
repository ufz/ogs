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
#include <set>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>

#include "IEquationIdStorage.h"

namespace DiscreteLib
{

/**
 * \brief Equation ID storage with random point or equation IDs
 *
 * This storage contains every pair of point ID and equation ID in a memory.
 * Hence it consumes more memory but is the most flexible to represent
 * arbitrary relationships between DoFs and equation ID, e.g. limited distribution
 * of some variables in a mesh.
 */
class RandomEquationIdStorage : public IEquationIdStorage
{
public:
    /**
     * constructor
     * @param sorted_vec_pt_id
     */
    explicit RandomEquationIdStorage(const std::vector<std::size_t> &sorted_vec_pt_id);

    /**
     * constructor
     * @param pt_id_start   Point ID start
     * @param n             Point count
     */
    RandomEquationIdStorage(std::size_t pt_id_start, std::size_t n);

    ///
    virtual ~RandomEquationIdStorage() {};

    ///
    virtual bool isSequential() const {return false;};

    /**
     *
     * @param pt_id
     * @return
     */
    virtual bool hasPoint(std::size_t pt_id) const;

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
     * @param pt_id
     * @param eqs_id
     */
    virtual void set(std::size_t pt_id, std::size_t eqs_id);

    /**
     *
     * @param dof_start
     * @param delta_per_pt
     * @return
     */
    virtual std::size_t setAll(std::size_t dof_start, std::size_t delta_per_pt);

    /**
     *
     * @return
     */
    virtual std::size_t size() const {return _map_pt2eqs.size();};

    /**
     *
     * @param pt_id
     * @return
     */
    virtual std::size_t equationID(std::size_t pt_id) const;

    /**
     *
     * @param address_id
     * @return
     */
    virtual std::size_t pointID(std::size_t address_id) const;

    /// print debug info
    virtual void printout() const;

private:
    typedef boost::bimaps::bimap<
                boost::bimaps::set_of<std::size_t>,
                boost::bimaps::multiset_of<std::size_t>
            > bimap_t;
    typedef bimap_t::value_type value_type;
    bimap_t _map_pt2eqs;
    std::vector<std::size_t> _list_pt_id;
    std::set<std::size_t> _deactive;
};

} //end


#endif //RANDOMEQUATIONIDSTORAGE_H_
