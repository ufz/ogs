/**
 * \file   IEquationIdStorage.h
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
#ifndef IEQUATIONIDSTORAGE_H_
#define IEQUATIONIDSTORAGE_H_

namespace DiscreteLib
{

/**
 * \brief Interface of a mapping table between point id and equation id
 */
class IEquationIdStorage
{
public:
    static const std::size_t index_npos = -1;

    /// destructor
    virtual ~IEquationIdStorage() {};

    /**
     * check if the given pointID is stored
     * @param pt_id Discrete point id
     * @return true/false
     */
    virtual bool hasPoint(std::size_t pt_id) const = 0;

    /**
     * check if the given value is stored
     * @param eqs_id    Equation id
     * @return  true/false
     */
    virtual bool hasEquationID(std::size_t eqs_id) const = 0;

    /**
     * get a range of keys
     * @param i_start   Start index
     * @param i_end     End index
     */
    virtual void getPointRange(std::size_t &i_start, std::size_t &i_end) const = 0;

    /**
     * activate a discrete point
     * @param pt_id     Discrete point id
     * @param b         True:active, False:inactive
     */
    virtual void activate(std::size_t pt_id, bool b) = 0;

    /**
     * return activeness of a given point
     * @param pt_id     Discrete point id
     * @return true/false
     */
    virtual bool isActive(std::size_t pt_id) const = 0;

    /**
     * set a mapping between a discrete point id and equation id
     * @param pt_id     Discrete point id
     * @param eqs_id    Equation id
     */
    virtual void set(std::size_t pt_id, std::size_t eqs_id) = 0;

    /**
     * renumber all equation indexes
     *
     * @param address_start     The beginning of the equation id
     * @param delta             Increment of equation id per point. Default is 1.
     * @return the last equation ID
     */
    virtual std::size_t setAll(std::size_t address_start, std::size_t delta=1) = 0;

    /**
     * get a size of stored points
     * @return
     */
    virtual std::size_t size() const = 0;

    /**
     * get an equationID of the given pointID
     * @param pt_id
     * @return
     *  index_npos (-1) is returned for invalid point id
     */
    virtual std::size_t equationID(std::size_t pt_id) const = 0;

    /**
     * get a pointID from the given equationID
     * @param eqs_id
     * @return
     *  index_npos (-1) is returned for invalid equation id
     */
    virtual std::size_t pointID(std::size_t eqs_id) const = 0;

    /**
     * return if this storage assumes sequential numbering of equation ID
     * @return
     */
    virtual bool isSequential() const = 0;

    /// print debug info
    virtual void printout() const = 0;
};

} //end

#endif //IEQUATIONIDSTORAGE_H_
