/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file OptionLeaf.h
 *
 * Created on 2012-02-02 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <ostream>
#include <sstream>
#include <limits>

#include "OptionNode.h"

namespace BaseLib
{

/**
 * \brief Leaf node in Options tree
 *
 * \tparam value type
 */
template <class T>
class OptionLeaf : public OptionNode
{
public:
    /**
     * Constructor
     * @param val   value
     */
    explicit OptionLeaf(const T &val) : _value(val) {};

    /**
     *
     */
    virtual ~OptionLeaf() {};

    /**
     * return if this node has value
     * @return true
     */
    bool isValue() const {return true;};

    /**
     * return if this node has string value
     * @return false (default)
     */
    virtual inline bool isString() const {return false;};

    /**
     * return value of this leaf
     * @return value
     */
    const T& getValue() const {return _value;};

    /**
     * return value as text data
     * @return string
     */
    virtual inline std::string getText() const
    {
        std::stringstream ss;
        ss << _value;
        return ss.str();
    };

    /**
     * print out content of this node
     *
     * @param os    output stream
     * @param depth depth of this node
     */
    virtual inline void printout(std::ostream &/*os*/, size_t /*depth*/) const {};

private:
    T _value;
};

template <>
inline bool OptionLeaf<std::string>::isString() const {return true;};

template <>
inline std::string OptionLeaf<std::string>::getText() const {return _value;};

template <>
inline void OptionLeaf<std::string>::printout(std::ostream &os, size_t depth) const
{
    for (size_t i=0; i<depth; i++)
        os << "\t";
    os << "value: " << _value << std::endl;
};

template <>
inline std::string OptionLeaf<double>::getText() const
{
    std::stringstream ss;
    ss.precision(std::numeric_limits<double>::digits10);
    ss << std::scientific << _value;
    return ss.str();
};

} // end BaseLib

