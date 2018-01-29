/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include <Eigen/Dense>

namespace MathLib
{
/// The invariants and the Kelving mapping are explained in detail in the
/// article "On Advantages of the Kelvin Mapping in Finite Element
/// Implementations of Deformation Processes" \cite Nagel2016.
namespace KelvinVector
{
/// Kelvin vector dimensions for given displacement dimension.
template <int DisplacementDim>
struct KelvinVectorDimensions;

template <>
struct KelvinVectorDimensions<2>
{
    static int const value = 4;
};

template <>
struct KelvinVectorDimensions<3>
{
    static int const value = 6;
};

//
// Kelvin vector and matrix templates for given displacement dimension.
//

/// Kelvin vector type for given displacement dimension.
/// \note The Eigen vector is always a fixed size vector in contrast to a shape
/// matrix policy types like BMatrixPolicyType::KelvinVectorType.
template <int DisplacementDim>
using KelvinVectorType =
    Eigen::Matrix<double,
                  KelvinVectorDimensions<DisplacementDim>::value,
                  1,
                  Eigen::ColMajor>;

/// Kelvin matrix type for given displacement dimension.
/// \note The Eigen matrix is always a fixed size vector in contrast to a shape
/// matrix policy types like BMatrixPolicyType::KelvinMatrixType.
template <int DisplacementDim>
using KelvinMatrixType =
    Eigen::Matrix<double,
                  KelvinVectorDimensions<DisplacementDim>::value,
                  KelvinVectorDimensions<DisplacementDim>::value,
                  Eigen::RowMajor>;

/// Invariants used in mechanics, based on Kelvin representation of the vectors
/// and matrices.
/// The invariants are computed at process creation time.
template <int KelvinVectorSize>
struct Invariants final
{
    static_assert(KelvinVectorSize == 4 || KelvinVectorSize == 6,
                  "KelvinVector invariants for vectors of size different than "
                  "4 or 6 is not allowed.");
    /// Kelvin mapping of deviatoric projection tensor. \f$A_{\rm dev} = P_{\rm
    /// dev}:A\f$ for \f$A\f$ being a second order tensor.
    static Eigen::Matrix<double, KelvinVectorSize, KelvinVectorSize> const
        deviatoric_projection;
    /// Kelvin mapping of spherical projection tensor. \f$A_{\rm sph} = P_{\rm
    /// sph}:A\f$ for \f$A\f$ being a second order tensor.
    static Eigen::Matrix<double, KelvinVectorSize, KelvinVectorSize> const
        spherical_projection;
    /// Kelvin mapping of 2nd order identity tensor.
    static Eigen::Matrix<double, KelvinVectorSize, 1> const identity2;

    /// Determinant of a matrix in Kelvin vector representation.
    static double determinant(
        Eigen::Matrix<double, KelvinVectorSize, 1> const& v);

    /// The von Mises equivalent stress.
    /// \note The input vector must have trace equal zero.
    static double equivalentStress(
        Eigen::Matrix<double, KelvinVectorSize, 1> const& deviatoric_v);

    /// Second invariant of deviatoric tensor.
    /// \note The input vector must have trace equal zero.
    static double J2(
        Eigen::Matrix<double, KelvinVectorSize, 1> const& deviatoric_v);

    /// Third invariant, equal to determinant of a deviatoric tensor.
    /// \note The input vector must have trace equal zero.
    static double J3(
        Eigen::Matrix<double, KelvinVectorSize, 1> const& deviatoric_v);

    /// Trace of the corresponding tensor.
    static double trace(Eigen::Matrix<double, KelvinVectorSize, 1> const& v);
};

//
// Inverses of a Kelvin vector.
//

/// Inverse of a matrix in Kelvin vector representation.
/// There are only implementations for the Kelvin vector size 4 and 6.
template <int KelvinVectorSize>
Eigen::Matrix<double, KelvinVectorSize, 1, Eigen::ColMajor, KelvinVectorSize, 1>
inverse(Eigen::Matrix<double,
                      KelvinVectorSize,
                      1,
                      Eigen::ColMajor,
                      KelvinVectorSize,
                      1> const& v);

/// Conversion of a Kelvin vector to a 3x3 matrix
/// Only implementations for KelvinVectorSize 4 and 6 are provided.
template <int KelvinVectorSize>
Eigen::Matrix<double, 3, 3> kelvinToTensor(
    Eigen::Matrix<double, KelvinVectorSize, 1, Eigen::ColMajor> const& v);

}  // namespace KelvinVector
}  // namespace MathLib

#include "KelvinVector-impl.h"
