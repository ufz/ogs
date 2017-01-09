/**
 * \file    ElementErrorCode.h
 * \author  Karsten Rink
 * \date    2014-02-21
 * \brief   Definition of ElementErrorCodes.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ELEMENTERRORCODES_H
#define ELEMENTERRORCODES_H

#include <bitset>
#include <string>


/// Possible error flags for mesh elements
enum class ElementErrorFlag
{
    ZeroVolume,
    NonCoplanar,
    NonConvex,
    NodeOrder,
    //... add other error flags here
    MaxValue // this needs to be last to set the bitset size correctly!
};

/// Collects error flags for mesh elements
class ElementErrorCode : public std::bitset<static_cast<std::size_t>(ElementErrorFlag::MaxValue)>
{
public:
    /// Get value for a specific flag
    bool get(ElementErrorFlag e) const { return test(static_cast<std::size_t>(e)); }
    /// Set a specific flag
    void set(ElementErrorFlag e) { std::bitset<static_cast<std::size_t>(ElementErrorFlag::MaxValue)>::set(static_cast<std::size_t>(e), true); }
    /// Reset a specific flag
    void reset(ElementErrorFlag e) { std::bitset<static_cast<std::size_t>(ElementErrorFlag::MaxValue)>::set(static_cast<std::size_t>(e), false); }

    inline reference operator[](const ElementErrorFlag e) { return std::bitset<static_cast<std::size_t>(ElementErrorFlag::MaxValue)>::operator[](static_cast<std::size_t>(e)); }
    inline bool operator[](const ElementErrorFlag e) const { return std::bitset<static_cast<std::size_t>(ElementErrorFlag::MaxValue)>::operator[](static_cast<std::size_t>(e)); }

    /// Returns a string output for a specific error flag
    static std::string toString(const ElementErrorFlag e)
    {
        if (e == ElementErrorFlag::ZeroVolume)
            return "zero volume";
        else if (e == ElementErrorFlag::NonCoplanar)
            return "non coplanar nodes";
        else if (e == ElementErrorFlag::NonConvex)
            return "non-convex geometry";
        else if (e == ElementErrorFlag::NodeOrder)
            return "wrong node order";
        return "nonspecified error";
    }

private:

};


#endif //ELEMENTERRORCODES_H
