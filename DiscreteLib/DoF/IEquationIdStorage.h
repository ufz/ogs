/**
 * \file   IEquationIdStorage.h
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

    virtual ~IEquationIdStorage() {};

    virtual bool hasKey(std::size_t pt_id) const = 0;

    virtual bool hasValue(std::size_t eqs_id) const = 0;

    virtual void key_range(std::size_t &i_start, std::size_t &i_end) const = 0;

    virtual void activate(std::size_t pt_id, bool b) = 0;

    virtual bool isActive(std::size_t pt_id) const = 0;

    virtual void set(std::size_t pt_id, long eqs_id) = 0;

    virtual std::size_t setAll(std::size_t address_start, std::size_t dn_pt=1) = 0;

    virtual std::size_t size() const = 0;

    virtual std::size_t address(std::size_t key_id) const = 0;

    virtual std::size_t key(std::size_t address_id) const = 0;

    virtual bool isSequential() const = 0;
};

} //end

#endif //IEQUATIONIDSTORAGE_H_
