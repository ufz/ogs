/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CoordinateSystem.h"

#include <typeinfo>

#include "Parameter.h"

namespace ProcessLib
{
CoordinateSystem::CoordinateSystem(Parameter<double> const& e0,
                                   Parameter<double> const& e1)
    : _base{&e0, &e1, nullptr}
{
    if (typeid(_base[0]) != typeid(_base[1]))
    {
        OGS_FATAL(
            "The parameter types for the basis must be equal but they are '%s' "
            "and '%s'.",
            typeid(_base[0]).name(),
            typeid(_base[1]).name());
    }
    if (_base[0]->isTimeDependent() || _base[1]->isTimeDependent())
    {
        OGS_FATAL("The parameters for the basis must not be time dependent.");
    }
    if (_base[0]->getNumberOfComponents() != 2 ||
        _base[1]->getNumberOfComponents() != 2)
    {
        OGS_FATAL("The parameters for the 2D basis must have two components.");
    }
}

CoordinateSystem::CoordinateSystem(Parameter<double> const& e0,
                                   Parameter<double> const& e1,
                                   Parameter<double> const& e2)
    : _base{&e0, &e1, &e2}
{
    if ((typeid(_base[0]) != typeid(_base[1])) ||
        (typeid(_base[1]) != typeid(_base[2])) ||
        (typeid(_base[2]) != typeid(_base[0])))
    {
        OGS_FATAL(
            "The parameter types for the basis must be equal but they are "
            "'%s', '%s', and '%s'.",
            typeid(_base[0]).name(),
            typeid(_base[1]).name(),
            typeid(_base[2]).name());
    }
    if (_base[0]->isTimeDependent() || _base[1]->isTimeDependent(),
        _base[2]->isTimeDependent())
    {
        OGS_FATAL("The parameters for the basis must not be time dependent.");
    }
    if (_base[0]->getNumberOfComponents() != 3 ||
        _base[1]->getNumberOfComponents() != 3 ||
        _base[2]->getNumberOfComponents() != 3)
    {
        OGS_FATAL(
            "The parameters for the 3D basis must have three components.");
    }
}

template <>
Eigen::Matrix<double, 2, 2> CoordinateSystem::transformation<2>(
    SpatialPosition const& pos) const
{
    if (_base[2] != nullptr)
    {
        OGS_FATAL(
            "The coordinate system is 3D but a transformation for 2D case is "
            "requested.");
    }

    auto e0 = (*_base[0])(0 /* time independent */, pos);
    auto e1 = (*_base[1])(0 /* time independent */, pos);
    Eigen::Matrix<double, 2, 2> t;
    t << e0[0], e1[0], e0[1], e1[1];

    return t;
}

template <>
Eigen::Matrix<double, 3, 3> CoordinateSystem::transformation<3>(
    SpatialPosition const& pos) const
{
    if (_base[2] == nullptr)
    {
        OGS_FATAL(
            "The coordinate system is 2D but a transformation for 3D case is "
            "requested.");
    }

    auto e0 = (*_base[0])(0 /* time independent */, pos);
    auto e1 = (*_base[1])(0 /* time independent */, pos);
    auto e2 = (*_base[2])(0 /* time independent */, pos);
    Eigen::Matrix<double, 3, 3> t;
    t << e0[0], e1[0], e2[0], e0[1], e1[1], e2[1], e0[2], e1[2], e2[2];

#ifndef NDEBUG
    if (std::abs(t.determinant() - 1) > std::numeric_limits<double>::epsilon())
    {
        OGS_FATAL(
            "The determinant of the coordinate system transformation matrix is "
            "'%g', which is not sufficiently close to unity.",
            t.determinant());
    }
#endif  // NDEBUG
    return t;
}
}  // namespace ProcessLib
