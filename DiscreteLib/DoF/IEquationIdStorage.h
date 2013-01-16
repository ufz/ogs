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
 * \brief Interface of mapped address
 */
class IEquationIdStorage
{
public:
    static const std::size_t index_npos = -1;

    /// destructor
    virtual ~IEquationIdStorage() {};

    /**
     * check if the given key is stored
     * @param pt_id Discrete point id
     * @return true/false
     */
    virtual bool hasKey(std::size_t pt_id) const = 0;

    /**
     * check if the given value is stored
     * @param eqs_id    Equation id
     * @return  true/false
     */
    virtual bool hasValue(std::size_t eqs_id) const = 0;

    /**
     * get a range of keys
     * @param i_start   Start index
     * @param i_end     End index
     */
    virtual void key_range(std::size_t &i_start, std::size_t &i_end) const = 0;

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
    virtual void set(std::size_t pt_id, long eqs_id) = 0;

    /**
     *
     * @param address_start
     * @param dn_pt
     * @return
     */
    virtual std::size_t setAll(std::size_t address_start, std::size_t dn_pt=1) = 0;

    /**
     * get a size of stored values
     * @return
     */
    virtual std::size_t size() const = 0;

    /**
     * get an address of the given key
     * @param key_id
     * @return
     */
    virtual std::size_t address(std::size_t key_id) const = 0;

    /**
     * get a key from the given address
     * @param address_id
     * @return
     */
    virtual std::size_t key(std::size_t address_id) const = 0;

    /**
     * return if this storage assumes sequential
     * @return
     */
    virtual bool isSequential() const = 0;
};

} //end

#endif //IEQUATIONIDSTORAGE_H_
