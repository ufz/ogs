/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "KelvinVector.h"

#include "BaseLib/Error.h"

namespace MathLib
{
namespace KelvinVector
{
template <>
double Invariants<6>::determinant(Eigen::Matrix<double, 6, 1> const& v)
{
    return v(0) * v(1) * v(2) + v(3) * v(4) * v(5) / std::sqrt(2.) -
           v(3) * v(3) * v(2) / 2. - v(4) * v(4) * v(0) / 2. -
           v(5) * v(5) * v(1) / 2.;
}

template <>
double Invariants<4>::determinant(Eigen::Matrix<double, 4, 1> const& v)
{
    return v(2) * (v(0) * v(1) - v(3) * v(3) / 2.);
}

template <>
Eigen::Matrix<double, 4, 1, Eigen::ColMajor, 4, 1> inverse(
    Eigen::Matrix<double, 4, 1, Eigen::ColMajor, 4, 1> const& v)
{
    assert(Invariants<4>::determinant(v) != 0);

    Eigen::Matrix<double, 4, 1, Eigen::ColMajor, 4, 1> inv;
    inv(0) = v(1) * v(2);
    inv(1) = v(0) * v(2);
    inv(2) = v(0) * v(1) - v(3) * v(3) / 2.;
    inv(3) = -v(3) * v(2);
    return inv / Invariants<4>::determinant(v);
}

template <>
Eigen::Matrix<double, 6, 1, Eigen::ColMajor, 6, 1> inverse(
    Eigen::Matrix<double, 6, 1, Eigen::ColMajor, 6, 1> const& v)
{
    assert(Invariants<6>::determinant(v) != 0);

    Eigen::Matrix<double, 6, 1, Eigen::ColMajor, 6, 1> inv;
    inv(0) = v(1) * v(2) - v(4) * v(4) / 2.;
    inv(1) = v(0) * v(2) - v(5) * v(5) / 2.;
    inv(2) = v(0) * v(1) - v(3) * v(3) / 2.;
    inv(3) = v(4) * v(5) / std::sqrt(2.) - v(3) * v(2);
    inv(4) = v(3) * v(5) / std::sqrt(2.) - v(4) * v(0);
    inv(5) = v(4) * v(3) / std::sqrt(2.) - v(1) * v(5);
    return inv / Invariants<6>::determinant(v);
}

template <>
Eigen::Matrix<double, 3, 3> kelvinVectorToTensor(
    Eigen::Matrix<double, 4, 1, Eigen::ColMajor, 4, 1> const& v)
{
    Eigen::Matrix<double, 3, 3> m;
    m << v[0], v[3] / std::sqrt(2.), 0, v[3] / std::sqrt(2.), v[1], 0, 0, 0,
        v[2];
    return m;
}

template <>
Eigen::Matrix<double, 3, 3> kelvinVectorToTensor(
    Eigen::Matrix<double, 6, 1, Eigen::ColMajor, 6, 1> const& v)
{
    Eigen::Matrix<double, 3, 3> m;
    m << v[0], v[3] / std::sqrt(2.), v[5] / std::sqrt(2.), v[3] / std::sqrt(2.),
        v[1], v[4] / std::sqrt(2.), v[5] / std::sqrt(2.), v[4] / std::sqrt(2.),
        v[2];
    return m;
}

template <>
Eigen::Matrix<double, 3, 3> kelvinVectorToTensor(Eigen::Matrix<double,
                                                               Eigen::Dynamic,
                                                               1,
                                                               Eigen::ColMajor,
                                                               Eigen::Dynamic,
                                                               1> const& v)
{
    if (v.size() == 4)
    {
        Eigen::Matrix<double, 4, 1, Eigen::ColMajor, 4, 1> v4;
        v4 << v[0], v[1], v[2], v[3];
        return kelvinVectorToTensor(v4);
    }
    if (v.size() == 6)
    {
        Eigen::Matrix<double, 6, 1, Eigen::ColMajor, 6, 1> v6;
        v6 << v[0], v[1], v[2], v[3], v[4], v[5];
        return kelvinVectorToTensor(v6);
    }
    OGS_FATAL(
        "Conversion of dynamic Kelvin vector of size {:d} to a tensor is not "
        "possible. Kelvin vector must be of size 4 or 6.",
        v.size());
}

template <>
KelvinVectorType<2> tensorToKelvin<2>(Eigen::Matrix<double, 3, 3> const& m)
{
    assert(std::abs(m(0, 1) - m(1, 0)) <
           std::numeric_limits<double>::epsilon());
    assert(m(0, 2) == 0);
    assert(m(1, 2) == 0);
    assert(m(2, 0) == 0);
    assert(m(2, 1) == 0);

    KelvinVectorType<2> v;
    v << m(0, 0), m(1, 1), m(2, 2), m(0, 1) * std::sqrt(2.);
    return v;
}

template <>
KelvinVectorType<3> tensorToKelvin<3>(Eigen::Matrix<double, 3, 3> const& m)
{
    assert(std::abs(m(0, 1) - m(1, 0)) <
           std::numeric_limits<double>::epsilon());
    assert(std::abs(m(1, 2) - m(2, 1)) <
           std::numeric_limits<double>::epsilon());
    assert(std::abs(m(0, 2) - m(2, 0)) <
           std::numeric_limits<double>::epsilon());

    KelvinVectorType<3> v;
    v << m(0, 0), m(1, 1), m(2, 2), m(0, 1) * std::sqrt(2.),
        m(1, 2) * std::sqrt(2.), m(0, 2) * std::sqrt(2.);
    return v;
}

template <>
Eigen::Matrix<double, 4, 1> kelvinVectorToSymmetricTensor(
    Eigen::Matrix<double, 4, 1, Eigen::ColMajor, 4, 1> const& v)
{
    Eigen::Matrix<double, 4, 1> m;
    m << v[0], v[1], v[2], v[3] / std::sqrt(2.);
    return m;
}

template <>
Eigen::Matrix<double, 6, 1> kelvinVectorToSymmetricTensor(
    Eigen::Matrix<double, 6, 1, Eigen::ColMajor, 6, 1> const& v)
{
    Eigen::Matrix<double, 6, 1> m;
    m << v[0], v[1], v[2], v[3] / std::sqrt(2.), v[4] / std::sqrt(2.),
        v[5] / std::sqrt(2.);
    return m;
}

template <>
Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor, Eigen::Dynamic, 1>
kelvinVectorToSymmetricTensor(Eigen::Matrix<double,
                                            Eigen::Dynamic,
                                            1,
                                            Eigen::ColMajor,
                                            Eigen::Dynamic,
                                            1> const& v)
{
    if (v.size() == 4)
    {
        return kelvinVectorToSymmetricTensor<4>(v);
    }
    if (v.size() == 6)
    {
        return kelvinVectorToSymmetricTensor<6>(v);
    }
    OGS_FATAL(
        "Kelvin vector to tensor conversion expected an input vector of size 4 "
        "or 6, but a vector of size {:d} was given.",
        v.size());
}

template <>
KelvinMatrixType<2> fourthOrderRotationMatrix<2>(
    Eigen::Matrix<double, 2, 2, Eigen::ColMajor, 2, 2> const& transformation)
{
    // 1-based index access for convenience.
    auto Q = [&](int const i, int const j) {
        return transformation(i - 1, j - 1);
    };

    MathLib::KelvinVector::KelvinMatrixType<2> R;
    R << Q(1, 1) * Q(1, 1), Q(1, 2) * Q(1, 2), 0,
        std::sqrt(2) * Q(1, 1) * Q(1, 2), Q(2, 1) * Q(2, 1), Q(2, 2) * Q(2, 2),
        0, std::sqrt(2) * Q(2, 1) * Q(2, 2), 0, 0, 1, 0,
        std::sqrt(2) * Q(1, 1) * Q(2, 1), std::sqrt(2) * Q(1, 2) * Q(2, 2), 0,
        Q(1, 1) * Q(2, 2) + Q(1, 2) * Q(2, 1);
    return R;
}

template <>
KelvinMatrixType<3> fourthOrderRotationMatrix<3>(
    Eigen::Matrix<double, 3, 3, Eigen::ColMajor, 3, 3> const& transformation)
{
    // 1-based index access for convenience.
    auto Q = [&](int const i, int const j) {
        return transformation(i - 1, j - 1);
    };

    MathLib::KelvinVector::KelvinMatrixType<3> R;
    R << Q(1, 1) * Q(1, 1), Q(1, 2) * Q(1, 2), Q(1, 3) * Q(1, 3),
        std::sqrt(2) * Q(1, 1) * Q(1, 2), std::sqrt(2) * Q(1, 2) * Q(1, 3),
        std::sqrt(2) * Q(1, 1) * Q(1, 3), Q(2, 1) * Q(2, 1), Q(2, 2) * Q(2, 2),
        Q(2, 3) * Q(2, 3), std::sqrt(2) * Q(2, 1) * Q(2, 2),
        std::sqrt(2) * Q(2, 2) * Q(2, 3), std::sqrt(2) * Q(2, 1) * Q(2, 3),
        Q(3, 1) * Q(3, 1), Q(3, 2) * Q(3, 2), Q(3, 3) * Q(3, 3),
        std::sqrt(2) * Q(3, 1) * Q(3, 2), std::sqrt(2) * Q(3, 2) * Q(3, 3),
        std::sqrt(2) * Q(3, 1) * Q(3, 3), std::sqrt(2) * Q(1, 1) * Q(2, 1),
        std::sqrt(2) * Q(1, 2) * Q(2, 2), std::sqrt(2) * Q(1, 3) * Q(2, 3),
        Q(1, 1) * Q(2, 2) + Q(1, 2) * Q(2, 1),
        Q(1, 2) * Q(2, 3) + Q(1, 3) * Q(2, 2),
        Q(1, 1) * Q(2, 3) + Q(1, 3) * Q(2, 1), std::sqrt(2) * Q(2, 1) * Q(3, 1),
        std::sqrt(2) * Q(2, 2) * Q(3, 2), std::sqrt(2) * Q(2, 3) * Q(3, 3),
        Q(2, 1) * Q(3, 2) + Q(2, 2) * Q(3, 1),
        Q(2, 2) * Q(3, 3) + Q(2, 3) * Q(3, 2),
        Q(2, 1) * Q(3, 3) + Q(2, 3) * Q(3, 1), std::sqrt(2) * Q(1, 1) * Q(3, 1),
        std::sqrt(2) * Q(1, 2) * Q(3, 2), std::sqrt(2) * Q(1, 3) * Q(3, 3),
        Q(1, 1) * Q(3, 2) + Q(1, 2) * Q(3, 1),
        Q(1, 2) * Q(3, 3) + Q(1, 3) * Q(3, 2),
        Q(1, 1) * Q(3, 3) + Q(1, 3) * Q(3, 1);
    return R;
}
}  // namespace KelvinVector
}  // namespace MathLib
