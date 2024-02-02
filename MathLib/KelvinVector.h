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

#include "BaseLib/Error.h"

namespace MathLib
{
/// The invariants and the Kelving mapping are explained in detail in the
/// article "On Advantages of the Kelvin Mapping in Finite Element
/// Implementations of Deformation Processes" \cite Nagel2016.
namespace KelvinVector
{
/// Kelvin vector dimensions for given displacement dimension.
constexpr int kelvin_vector_dimensions(int const displacement_dim)
{
    if (displacement_dim == 2)
    {
        return 4;
    }
    else if (displacement_dim == 3)
    {
        return 6;
    }
    OGS_FATAL(
        "Cannot convert displacement dimension {} to kelvin vector dimension.",
        displacement_dim);
}

//
// Kelvin vector and matrix templates for given displacement dimension.
//

/// Kelvin vector type for given displacement dimension.
/// \note The Eigen vector is always a fixed size vector in contrast to a shape
/// matrix policy types like BMatrixPolicyType::KelvinVectorType.
template <int DisplacementDim>
using KelvinVectorType =
    Eigen::Matrix<double, kelvin_vector_dimensions(DisplacementDim), 1,
                  Eigen::ColMajor>;

/// Kelvin matrix type for given displacement dimension.
/// \note The Eigen matrix is always a fixed size vector in contrast to a shape
/// matrix policy types like BMatrixPolicyType::KelvinMatrixType.
template <int DisplacementDim>
using KelvinMatrixType =
    Eigen::Matrix<double, kelvin_vector_dimensions(DisplacementDim),
                  kelvin_vector_dimensions(DisplacementDim), Eigen::RowMajor>;

/// Returns an expressions for a Kelvin vector filled with zero.
template <int DisplacementDim>
constexpr auto KVzero()
{
    return KelvinVectorType<DisplacementDim>::Zero();
}

/// Returns an expressions for a Kelvin matrix filled with zero.
template <int DisplacementDim>
constexpr auto KMzero()
{
    return KelvinMatrixType<DisplacementDim>::Zero();
}

/// Returns an expressions for a Kelvin vector filled with NaN.
template <int DisplacementDim>
constexpr auto KVnan()
{
    return KelvinVectorType<DisplacementDim>::Constant(
        std::numeric_limits<double>::quiet_NaN());
}

/// Returns an expressions for a Kelvin matrix filled with NaN.
template <int DisplacementDim>
constexpr auto KMnan()
{
    return KelvinMatrixType<DisplacementDim>::Constant(
        std::numeric_limits<double>::quiet_NaN());
}

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

    /// Get the norm of the deviatoric stress.
    static double FrobeniusNorm(
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

    /// Diagonal of the corresponding tensor which is always of length 3 in 2D
    /// and 3D cases.
    static Eigen::Vector3d diagonal(
        Eigen::Matrix<double, KelvinVectorSize, 1> const& v);
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
Eigen::Matrix<double, 3, 3> kelvinVectorToTensor(Eigen::Matrix<double,
                                                               KelvinVectorSize,
                                                               1,
                                                               Eigen::ColMajor,
                                                               KelvinVectorSize,
                                                               1> const& v);

/// Conversion of a 3x3 matrix to a Kelvin vector.
/// Only implementations for KelvinVectorSize 4 and 6 are provided.
template <int DisplacementDim>
KelvinVectorType<DisplacementDim> tensorToKelvin(
    Eigen::Matrix<double, 3, 3> const& m);

/// Conversion of a Kelvin vector to a short vector representation of a
/// symmetric 3x3 matrix.
///
/// In the 2D case the entries for the xx, yy, zz, and xy components are stored.
/// In the 3D case the entries for the xx, yy, zz, xy, yz, and xz components in
/// that particular order are stored.
///
/// This is opposite of the symmetricTensorToKelvinVector()
///
/// Only implementations for KelvinVectorSize 4 and 6, and dynamic size vectors
/// are provided.
template <int KelvinVectorSize>
Eigen::Matrix<double, KelvinVectorSize, 1, Eigen::ColMajor, KelvinVectorSize, 1>
kelvinVectorToSymmetricTensor(Eigen::Matrix<double,
                                            KelvinVectorSize,
                                            1,
                                            Eigen::ColMajor,
                                            KelvinVectorSize,
                                            1> const& v);

/// Conversion of a short vector representation of a
/// symmetric 3x3 matrix to a Kelvin vector.
///
/// This is opposite of the kelvinVectorToSymmetricTensor()
///
/// Only implementations for KelvinVectorSize 4 and 6, and dynamic size vectors
/// are provided.
template <typename Derived>
Eigen::Matrix<double, Eigen::MatrixBase<Derived>::RowsAtCompileTime, 1>
symmetricTensorToKelvinVector(Eigen::MatrixBase<Derived> const& v)
{
    static_assert(
        (Eigen::MatrixBase<Derived>::ColsAtCompileTime == 1) ||
            (Eigen::MatrixBase<Derived>::ColsAtCompileTime == Eigen::Dynamic),
        "KelvinVector must be a column vector");
    if (v.cols() != 1)
    {
        OGS_FATAL(
            "KelvinVector must be a column vector, but input has {:d} columns.",
            v.cols());
    }

    Eigen::Matrix<double, Eigen::MatrixBase<Derived>::RowsAtCompileTime, 1>
        result;
    if (v.rows() == 4)
    {
        result.resize(4, 1);
        result << v[0], v[1], v[2], v[3] * std::sqrt(2.);
    }
    else if (v.rows() == 6)
    {
        result.resize(6, 1);
        result << v[0], v[1], v[2], v[3] * std::sqrt(2.), v[4] * std::sqrt(2.),
            v[5] * std::sqrt(2.);
    }
    else
    {
        OGS_FATAL(
            "Symmetric tensor to Kelvin vector conversion expected an input "
            "vector of size 4 or 6, but a vector of size {:d} was given.",
            v.size());
    }
    return result;
}

/// Conversion of a short vector representation of a symmetric 3x3 matrix to a
/// Kelvin vector.
///
/// This overload takes a std::vector for the tensor values.
template <int DisplacementDim>
KelvinVectorType<DisplacementDim> symmetricTensorToKelvinVector(
    std::vector<double> const& values)
{
    constexpr int kelvin_vector_size =
        kelvin_vector_dimensions(DisplacementDim);

    if (values.size() != kelvin_vector_size)
    {
        OGS_FATAL(
            "Symmetric tensor to Kelvin vector conversion expected an input "
            "vector of size {:d}, but a vector of size {:d} was given.",
            kelvin_vector_size, values.size());
    }

    return symmetricTensorToKelvinVector(
        Eigen::Map<typename MathLib::KelvinVector::KelvinVectorType<
            DisplacementDim> const>(
            values.data(), kelvin_vector_dimensions(DisplacementDim), 1));
}

/// Lifting of a vector to a Kelvin matrix.
/// \f$ a -> A\f$ s.t. \f$k_{ij} a_{j} = A_[\alpha\beta} k_{\beta} \f$.
/// Conversion for 2D -> 4D and 3D -> 6D are implemented.
template <int DisplacementDim>
Eigen::Matrix<double, DisplacementDim,
              kelvin_vector_dimensions(DisplacementDim)>
liftVectorToKelvin(Eigen::Matrix<double, DisplacementDim, 1> const& v);

/// Reducing a Kelvin matrix to a vector.
/// Conversion for 4D -> 2D and 6D -> 3D are implemented.
template <int DisplacementDim>
Eigen::Matrix<double, DisplacementDim, 1> reduceKelvinToVector(
    Eigen::Matrix<double, DisplacementDim,
                  kelvin_vector_dimensions(DisplacementDim)> const& m);

/// Rotation tensor for Kelvin mapped vectors and tensors. It is meant to be
/// used for rotation of stress/strain tensors epsilon:Q and tangent stiffness
/// tensors Q*C*Q^t.
/// 2D and 3D implementations available.
template <int DisplacementDim>
KelvinMatrixType<DisplacementDim> fourthOrderRotationMatrix(
    Eigen::Matrix<double, DisplacementDim, DisplacementDim, Eigen::ColMajor,
                  DisplacementDim, DisplacementDim> const& transformation);

}  // namespace KelvinVector
}  // namespace MathLib

#include "KelvinVector-impl.h"
