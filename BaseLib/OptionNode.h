/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file OptionNode.h
 *
 * Created on 2012-02-02 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <ostream>

namespace BaseLib
{

/**
 * \brief Abstract class of nodes in Options tree
 */
class OptionNode 
{
public:
    /**
     *
     */
    OptionNode() {};

    /**
     *
     */
    virtual ~OptionNode() {};

    /**
     * return if this node has value
     * @return true if this node has value
     */
    virtual bool isValue() const {return false;};

    /**
     * return if this node has string value
     * @return false (default)
     */
    virtual bool isString() const {return false;};

    /**
     * return characters
     *
     * non-character data should be converted to string. Node without values
     * return empty string.
     * @return string
     */
    virtual std::string getText() const {return "";};

    /**
     * print out content
     *
     * @param os    output stream
     * @param depth depth of this node from certain position
     */
    virtual void printout(std::ostream &os, size_t depth=0) const = 0;

protected:
    /**
     * return dummy value in case there is no value
     *
     * @return
     */
    template<typename T>
    inline T getDummy() const {return 0;};

private:
    const std::string _dummy;
};

template<>
inline std::string OptionNode::getDummy<std::string>() const
{
    return _dummy;
};

} // end BaseLib
