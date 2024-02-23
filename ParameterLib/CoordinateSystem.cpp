/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CoordinateSystem.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <algorithm>
#include <iterator>
#include <limits>
#include <typeinfo>

#include "Parameter.h"

namespace ParameterLib
{
static double const tolerance = std::numeric_limits<double>::epsilon();
#ifndef NDEBUG
static constexpr char error_info[] =
    "The determinant of the coordinate system transformation matrix is '{:g}', "
    "which is not sufficiently close to unity with the tolerance of '{:g}'. "
    "Please adjust the accuracy of the local system bases";

static constexpr char normalization_error_info[] =
    "The direction vector given by parameter {:s} for local_coordinate_system "
    "is not normalized to unit length";
#endif  // NDEBUG

CoordinateSystem::CoordinateSystem(Parameter<double> const& unit_direction)
    : _base{nullptr, nullptr, &unit_direction}, _has_implicit_base(true)
{
    if (_base[2]->isTimeDependent())
    {
        OGS_FATAL(
            "The unit_normal parameter named {} must not be time dependent.",
            unit_direction.name);
    }
}

CoordinateSystem::CoordinateSystem(Parameter<double> const& e0,
                                   Parameter<double> const& e1)
    : _base{&e0, &e1, nullptr}, _has_implicit_base(false)
{
    if (typeid(_base[0]) != typeid(_base[1]))
    {
        OGS_FATAL(
            "The parameter types for the basis must be equal but they are "
            "'{:s}' and '{:s}'.",
            typeid(_base[0]).name(),
            typeid(_base[1]).name());
    }
    if (_base[0]->isTimeDependent() || _base[1]->isTimeDependent())
    {
        OGS_FATAL("The parameters for the basis must not be time dependent.");
    }
    if (_base[0]->getNumberOfGlobalComponents() != 2 ||
        _base[1]->getNumberOfGlobalComponents() != 2)
    {
        OGS_FATAL("The parameters for the 2D basis must have two components.");
    }
}

CoordinateSystem::CoordinateSystem(Parameter<double> const& e0,
                                   Parameter<double> const& e1,
                                   Parameter<double> const& e2)
    : _base{&e0, &e1, &e2}, _has_implicit_base(false)
{
    if ((typeid(_base[0]) != typeid(_base[1])) ||
        (typeid(_base[1]) != typeid(_base[2])) ||
        (typeid(_base[2]) != typeid(_base[0])))
    {
        OGS_FATAL(
            "The parameter types for the basis must be equal but they are "
            "'{:s}', '{:s}', and '{:s}'.",
            typeid(_base[0]).name(),
            typeid(_base[1]).name(),
            typeid(_base[2]).name());
    }
    if (_base[0]->isTimeDependent() || _base[1]->isTimeDependent(),
        _base[2]->isTimeDependent())
    {
        OGS_FATAL("The parameters for the basis must not be time dependent.");
    }
    if (_base[0]->getNumberOfGlobalComponents() != 3 ||
        _base[1]->getNumberOfGlobalComponents() != 3 ||
        _base[2]->getNumberOfGlobalComponents() != 3)
    {
        OGS_FATAL(
            "The parameters for the 3D basis must have three components.");
    }
}

Eigen::Matrix<double, 2, 2> getTransformationFromSingleBase2D(
    Parameter<double> const& unit_direction, SpatialPosition const& pos)
{
    Eigen::Matrix<double, 2, 2> t;
#ifndef NDEBUG
    if (std::abs(Eigen::Map<Eigen::Vector2d>(
                     unit_direction(0 /* time independent */, pos).data())
                     .norm() -
                 1) > tolerance)
    {
        OGS_FATAL(normalization_error_info, unit_direction.name);
    }
#endif  // NDEBUG

    auto normal = unit_direction(0 /* time independent */, pos);
    // base 0: ( normal[1], -normal[0])^T
    // base 1: ( normal[0], normal[1])^T
    t << normal[1], normal[0], -normal[0], normal[1];

#ifndef NDEBUG
    if (std::abs(t.determinant() - 1) > tolerance)
    {
        OGS_FATAL(error_info, t.determinant(), tolerance);
    }
#endif  // NDEBUG
    return t;
}

template <>
Eigen::Matrix<double, 2, 2> CoordinateSystem::transformation<2>(
    SpatialPosition const& pos) const
{
    if (_has_implicit_base)
    {
        return getTransformationFromSingleBase2D(*_base[2], pos);
    }

    if (_base[2] != nullptr)
    {
        OGS_FATAL(
            "The coordinate system is 3D but a transformation for 2D case is "
            "requested.");
    }

    Eigen::Matrix<double, 2, 2> t;
    t.col(0) = Eigen::Map<Eigen::Vector2d>(
        (*_base[0])(0 /* time independent */, pos).data());
    t.col(1) = Eigen::Map<Eigen::Vector2d>(
        (*_base[1])(0 /* time independent */, pos).data());

#ifndef NDEBUG
    if (std::abs(t.determinant() - 1) > tolerance)
    {
        OGS_FATAL(error_info, t.determinant(), tolerance);
    }
#endif  // NDEBUG
    return t;
}

Eigen::Matrix<double, 3, 3> getTransformationFromSingleBase3D(
    Parameter<double> const& unit_direction, SpatialPosition const& pos)
{
    auto e2 = unit_direction(0 /* time independent */, pos);
#ifndef NDEBUG
    if (std::abs(Eigen::Map<Eigen::Vector3d>(e2.data()).norm() - 1.0) >
        tolerance)
    {
        OGS_FATAL(normalization_error_info, unit_direction.name);
    }
#endif

    // Find the id of the first non-zero component of e2:
    auto const it = std::max_element(e2.begin(), e2.end(),
                                     [](const double a, const double b)
                                     { return std::abs(a) < std::abs(b); });
    const auto id = std::distance(e2.begin(), it);
    // Get other two component ids:
    const auto id_a = (id + 1) % 3;
    const auto id_b = (id + 2) % 3;

    // Compute basis vector e1 orthogonal to e2
    using Vector3 = Eigen::Vector3d;
    Vector3 e1(0.0, 0.0, 0.0);

    if (std::abs(e2[id_a]) < tolerance)
    {
        e1[id_a] = 1.0;
    }
    else if (std::abs(e2[id_b]) < tolerance)
    {
        e1[id_b] = 1.0;
    }
    else
    {
        e1[id_a] = 1.0;
        e1[id_b] = -e2[id_a] / e2[id_b];
    }

    e1.normalize();

    auto eigen_mapped_e2 = Eigen::Map<Vector3>(e2.data());

    auto e0 = e1.cross(eigen_mapped_e2);
    // |e0| = |e1 x e2| = |e1||e2|sin(theta) with theta the angle between e1 and
    // e2. Since |e1| = |e2| = 1.0, and theta = pi/2, we have |e0|=1. Therefore
    // e0 is normalized by nature.

    Eigen::Matrix<double, 3, 3> t;
    t.col(0) = e0;
    t.col(1) = e1;
    t.col(2) = eigen_mapped_e2;

#ifndef NDEBUG

    if (std::abs(t.determinant() - 1) > tolerance)
    {
        OGS_FATAL(error_info, t.determinant(), tolerance);
    }

#endif  // NDEBUG

    return t;
}

template <>
Eigen::Matrix<double, 3, 3> CoordinateSystem::transformation<3>(
    SpatialPosition const& pos) const
{
    if (_has_implicit_base)
    {
        return getTransformationFromSingleBase3D(*_base[2], pos);
    }

    if (_base[2] == nullptr)
    {
        OGS_FATAL(
            "The coordinate system is 2D but a transformation for 3D case is "
            "requested.");
    }

    Eigen::Matrix<double, 3, 3> t;
    t.col(0) = Eigen::Map<Eigen::Vector3d>(
        (*_base[0])(0 /* time independent */, pos).data());
    t.col(1) = Eigen::Map<Eigen::Vector3d>(
        (*_base[1])(0 /* time independent */, pos).data());
    t.col(2) = Eigen::Map<Eigen::Vector3d>(
        (*_base[2])(0 /* time independent */, pos).data());

#ifndef NDEBUG
    if (std::abs(t.determinant() - 1) > tolerance)
    {
        OGS_FATAL(error_info, t.determinant(), tolerance);
    }
#endif  // NDEBUG
    return t;
}

Eigen::Matrix<double, 3, 3> CoordinateSystem::transformationFromSingleBase_3d(
    SpatialPosition const& pos) const
{
    // If base 2, which stores the unit direction vector, has three components:
    if ((*_base[2])(0 /* time independent */, pos).size() == 3)
    {
        return getTransformationFromSingleBase3D(*_base[2], pos);
    }

    // If base 2, which stores the unit direction vector, has two components
    Eigen::Matrix<double, 3, 3> t = Eigen::Matrix<double, 3, 3>::Identity();
    t.template topLeftCorner<2, 2>() = transformation<2>(pos);

    return t;
}

Eigen::Matrix<double, 3, 3> CoordinateSystem::transformation_3d(
    SpatialPosition const& pos) const
{
    if (_has_implicit_base)
    {
        return transformationFromSingleBase_3d(pos);
    }

    if (_base[2] != nullptr)
    {
        return transformation<3>(pos);
    }

    auto e0 = (*_base[0])(0 /* time independent */, pos);
    auto e1 = (*_base[1])(0 /* time independent */, pos);
    Eigen::Matrix<double, 3, 3> t = Eigen::Matrix<double, 3, 3>::Identity();
    t.template topLeftCorner<2, 2>() << e0[0], e1[0], e0[1], e1[1];

#ifndef NDEBUG
    if (std::abs(t.determinant() - 1) > tolerance)
    {
        OGS_FATAL(error_info, t.determinant(), tolerance);
    }
#endif  // NDEBUG
    return t;
}

template <int Dimension>
Eigen::Matrix<double, Dimension, Dimension> CoordinateSystem::rotateTensor(
    std::vector<double> const& values, SpatialPosition const& pos) const
{
    assert(values.size() == Dimension * Dimension ||
           "Input vector has wrong dimension; expected 4 or 9 entries.");
    auto const tensor =
        Eigen::Map<Eigen::Matrix<double, Dimension, Dimension> const>(
            values.data(), Dimension, Dimension);
    auto const R = transformation<Dimension>(pos);
    return R * tensor * R.transpose();
}

template <int Dimension>
Eigen::Matrix<double, Dimension, Dimension>
CoordinateSystem::rotateDiagonalTensor(std::vector<double> const& values,
                                       SpatialPosition const& pos) const
{
    assert(values.size() == Dimension ||
           "Input vector has wrong dimension; expected 2 or 3 entries.");
    auto const tensor = Eigen::Map<Eigen::Matrix<double, Dimension, 1> const>(
        values.data(), Dimension, 1);
    auto const R = transformation<Dimension>(pos);
    return R * tensor.asDiagonal() * R.transpose();
}

template Eigen::Matrix<double, 2, 2> CoordinateSystem::rotateTensor<2>(
    std::vector<double> const& values, SpatialPosition const& pos) const;
template Eigen::Matrix<double, 3, 3> CoordinateSystem::rotateTensor<3>(
    std::vector<double> const& values, SpatialPosition const& pos) const;
template Eigen::Matrix<double, 2, 2> CoordinateSystem::rotateDiagonalTensor<2>(
    std::vector<double> const& values, SpatialPosition const& pos) const;
template Eigen::Matrix<double, 3, 3> CoordinateSystem::rotateDiagonalTensor<3>(
    std::vector<double> const& values, SpatialPosition const& pos) const;
}  // namespace ParameterLib
