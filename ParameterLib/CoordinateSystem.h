/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Core>
#include <array>
#include <vector>

namespace ParameterLib
{
template <typename T>
struct Parameter;
class SpatialPosition;
}  // namespace ParameterLib

namespace ParameterLib
{
/**
 * \brief A local coordinate system used for tensor transformations.
 *
 * It offers a simple way for input of anisotropic tensors w.r.t. a coordinate
 * system.
 * The basis vectors form a transformation matrix \f$R = (e_0, e_1, e_2)\f$.
 * For a given anisotropic tensor \f$A\f$ parameter with the corresponding
 * \ref ogs_file_param__prj__parameters__parameter__use_local_coordinate_system
 * the tensor is rotated according to the formula: \f$A' = R\cdot A\cdot R^T\f$.
 *
 * For computations in transverse isotropic material models, we can create a
 * coordinate system with only one base, where the last base is explicitly
 * given. The other bases are set as implicit and computed from the given base
 * as follows:
 *  - For a 2D coordinate system, the unit vector orthogonal to
 *    the given base is used as the first base,
 *  - For a 3D coordinate system, the given base vector, \c unit_direction,
 *    is set as the third base,  \f${\vec e}_2\f$. An arbitrary unit vector
 *    orthogonal to \f${\vec e}_2\f$ is selected as the second base
 *    \f$e_1\f$, and the first  base \f${\vec e}_0\f$ is calculated as
 *    \f${\vec e}_0 = {\vec e}_1 \times {\vec e}_2\f$.
 */
struct CoordinateSystem final
{
    /**
     * It is used to create a local coordinate system with only one base, where
     * the last base is explicitly given as \c unit_direction.
     *
     * @param unit_direction The specified unit direction.
     */
    explicit CoordinateSystem(Parameter<double> const& unit_direction);

    CoordinateSystem(Parameter<double> const& e0, Parameter<double> const& e1);

    CoordinateSystem(Parameter<double> const& e0,
                     Parameter<double> const& e1,
                     Parameter<double> const& e2);

    template <int Dimension>
    Eigen::Matrix<double, Dimension, Dimension> transformation(
        SpatialPosition const& pos) const;

    Eigen::Matrix<double, 3, 3> transformation_3d(
        SpatialPosition const& pos) const;

    template <int Dimension>
    Eigen::Matrix<double, Dimension, Dimension> rotateTensor(
        std::vector<double> const& values, SpatialPosition const& pos) const;

    template <int Dimension>
    Eigen::Matrix<double, Dimension, Dimension> rotateDiagonalTensor(
        std::vector<double> const& values, SpatialPosition const& pos) const;

private:
    std::array<Parameter<double> const*, 3> _base;
    bool _has_implicit_base;

    Eigen::Matrix<double, 3, 3> transformationFromSingleBase_3d(
        SpatialPosition const& pos) const;
};

extern template Eigen::Matrix<double, 2, 2> CoordinateSystem::rotateTensor<2>(
    std::vector<double> const& values, SpatialPosition const& pos) const;
extern template Eigen::Matrix<double, 3, 3> CoordinateSystem::rotateTensor<3>(
    std::vector<double> const& values, SpatialPosition const& pos) const;
extern template Eigen::Matrix<double, 2, 2>
CoordinateSystem::rotateDiagonalTensor<2>(std::vector<double> const& values,
                                          SpatialPosition const& pos) const;
extern template Eigen::Matrix<double, 3, 3>
CoordinateSystem::rotateDiagonalTensor<3>(std::vector<double> const& values,
                                          SpatialPosition const& pos) const;
}  // namespace ParameterLib
