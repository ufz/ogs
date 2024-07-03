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
#include <cmath>
#include <limits>
#include <typeinfo>

#include "MathLib/FormattingUtils.h"
#include "Parameter.h"

namespace ParameterLib
{
static double const tolerance = std::numeric_limits<double>::epsilon();

template <int Dim>
static void checkTransformationIsSON(
    Eigen::Matrix<double, Dim, Dim, Eigen::ColMajor, Dim, Dim> const& t)
{
    if (std::abs(t.determinant() - 1) > tolerance)
    {
        OGS_FATAL(
            "The determinant of the coordinate system transformation matrix is "
            "'{:g}', which is not sufficiently close to unity with the "
            "tolerance of '{:g}'. Please adjust the accuracy of the local "
            "system bases",
            t.determinant(), tolerance);
    }
    if (((t * t.transpose() - Eigen::Matrix<double, Dim, Dim>::Identity())
             .array() > 2 * tolerance)
            .any())
    {
        OGS_FATAL(
            "The transformation is not orthogonal because the difference "
            "A*A^T - I is:\n{}\nand at least one component deviates from zero "
            "by more then '{:g}'.",
            (t * t.transpose() - Eigen::Matrix<double, Dim, Dim>::Identity())
                .eval(),
            2 * tolerance);
    }
}

template <typename Derived>
static void checkNormalization(Eigen::MatrixBase<Derived> const& vec,
                               std::string_view const parmeter_name)
{
    if (std::abs(vec.squaredNorm() - 1.0) > tolerance)
    {
        OGS_FATAL(
            "The direction vector given by parameter {:s} for the "
            "local_coordinate_system is not normalized to unit length. Vector "
            "norm is {:g} and |v|^2-1 = {:g}.",
            parmeter_name, vec.norm(), vec.squaredNorm() - 1.0);
    }
}

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
    auto const& normal = unit_direction(0 /* time independent */, pos);
    checkNormalization(Eigen::Map<Eigen::Vector2d const>(normal.data()),
                       unit_direction.name);

    Eigen::Matrix<double, 2, 2> t;
    // base 0: ( normal[1], -normal[0])^T
    // base 1: ( normal[0], normal[1])^T
    t << normal[1], normal[0], -normal[0], normal[1];

    checkTransformationIsSON(t);
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

    checkTransformationIsSON(t);
    return t;
}

Eigen::Matrix<double, 3, 3> getTransformationFromSingleBase3D(
    Parameter<double> const& unit_direction, SpatialPosition const& pos)
{
    auto const& normal = unit_direction(0 /* time independent */, pos);

    Eigen::Matrix<double, 3, 3> t;
    auto e2 = t.col(2);
    e2 = Eigen::Map<Eigen::Vector3d const>(normal.data());
    checkNormalization(e2, unit_direction.name);

    // Find the id of the first non-zero component of e2:
    int id;
    e2.cwiseAbs().maxCoeff(&id);

    // Get other two component ids:
    const auto id_a = (id + 1) % 3;
    const auto id_b = (id + 2) % 3;

    // Compute basis vector e1 orthogonal to e2
    auto e1 = t.col(1);
    e1 = Eigen::Vector3d::Zero();

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

    // |e0| = |e1 x e2| = |e1||e2|sin(theta) with theta the angle between e1 and
    // e2. Since |e1| = |e2| = 1.0, and theta = pi/2, we have |e0|=1. Therefore
    // e0 is normalized by nature.
    t.col(0) = e1.cross(e2);

    checkTransformationIsSON(t);

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

    checkTransformationIsSON(t);
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

    checkTransformationIsSON(t);
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
    return rotateTensor<Dimension>(tensor, pos);
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
