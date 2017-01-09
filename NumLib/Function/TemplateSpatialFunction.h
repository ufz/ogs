/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TEMPLATESPATIALFUNCTION_H_
#define TEMPLATESPATIALFUNCTION_H_

#include "ISpatialFunction.h"

namespace NumLib
{

/**
 * Template class to relate pure mathematical functions with ISpatialFunction
 *
 * @tparam T_FUNCTION Function object which takes "double const* const" as an argument.
 */
template <class T_FUNCTION>
class TemplateSpatialFunction : public ISpatialFunction
{
public:
    /**
     * Constructor
     * @param f  a function object
     */
    TemplateSpatialFunction(const T_FUNCTION &f)
    : _f(f) {}

    /**
     * evaluate a function
     * @param pnt  a point object
     * @return evaluated value
     */
    virtual double operator()(const MathLib::Point3d& pnt) const
    {
        return _f(pnt.getCoords());
    }

private:
    /// a function object
    const T_FUNCTION _f;
};

}

#endif /* TEMPLATESPATIALFUNCTION_H_ */
