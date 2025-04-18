/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cmath>

namespace MathLib
{
namespace KelvinVector
{
template <int KelvinVectorSize>
double Invariants<KelvinVectorSize>::equivalentStress(
    Eigen::Matrix<double, KelvinVectorSize, 1> const& deviatoric_v)
{
    assert(std::abs(trace(deviatoric_v)) <=
           4e-12 * diagonal(deviatoric_v).norm());
    return std::sqrt(3 * J2(deviatoric_v));
}

template <int KelvinVectorSize>
double Invariants<KelvinVectorSize>::FrobeniusNorm(
    Eigen::Matrix<double, KelvinVectorSize, 1> const& deviatoric_v)
{
    return std::sqrt(deviatoric_v.transpose() * deviatoric_v);
}

template <int KelvinVectorSize>
double Invariants<KelvinVectorSize>::J2(
    Eigen::Matrix<double, KelvinVectorSize, 1> const& deviatoric_v)
{
    assert(std::abs(trace(deviatoric_v)) <=
           4e-12 * diagonal(deviatoric_v).norm());
    return 0.5 * deviatoric_v.transpose() * deviatoric_v;
}

/// Third invariant, equal to determinant of a deviatoric tensor.
/// \note The input vector must have trace equal zero.
template <int KelvinVectorSize>
double Invariants<KelvinVectorSize>::J3(
    Eigen::Matrix<double, KelvinVectorSize, 1> const& deviatoric_v)
{
    assert(std::abs(trace(deviatoric_v)) <=
           4e-12 * diagonal(deviatoric_v).norm());
    return determinant(deviatoric_v);
}

template <int KelvinVectorSize>
Eigen::Vector3d Invariants<KelvinVectorSize>::diagonal(
    Eigen::Matrix<double, KelvinVectorSize, 1> const& v)
{
    return v.template topLeftCorner<3, 1>();
}

/// Trace of the corresponding tensor.
template <int KelvinVectorSize>
double Invariants<KelvinVectorSize>::trace(
    Eigen::Matrix<double, KelvinVectorSize, 1> const& v)
{
    return diagonal(v).sum();
}

//
// Initialization of static Invariant variables.
//

namespace KelvinVector_detail
{
template <int KelvinVectorSize>
Eigen::Matrix<double, KelvinVectorSize, KelvinVectorSize>
initDeviatoricProjection()
{
    Eigen::Matrix<double, KelvinVectorSize, KelvinVectorSize> P_dev =
        Eigen::Matrix<double, KelvinVectorSize, KelvinVectorSize>::Identity();

    P_dev.template topLeftCorner<3, 3>() -=
        1. / 3. * Eigen::Matrix<double, 3, 3>::Ones();
    return P_dev;
}

template <int KelvinVectorSize>
Eigen::Matrix<double, KelvinVectorSize, KelvinVectorSize>
initSphericalProjection()
{
    Eigen::Matrix<double, KelvinVectorSize, KelvinVectorSize> P_sph =
        Eigen::Matrix<double, KelvinVectorSize, KelvinVectorSize>::Zero();

    P_sph.template topLeftCorner<3, 3>().setConstant(1. / 3.);
    return P_sph;
}

template <int KelvinVectorSize>
Eigen::Matrix<double, KelvinVectorSize, 1> initIdentity2()
{
    Eigen::Matrix<double, KelvinVectorSize, 1> ivec =
        Eigen::Matrix<double, KelvinVectorSize, 1>::Zero();

    ivec.template topLeftCorner<3, 1>().setConstant(1.);
    return ivec;
}

template <int KelvinVectorSize>
Eigen::Matrix<double, KelvinVectorSize, 1> initOnes2()
{
    Eigen::Matrix<double, KelvinVectorSize, 1> ivec =
        Eigen::Matrix<double, KelvinVectorSize, 1>::Ones();

    ivec.template bottomLeftCorner<KelvinVectorSize - 3, 1>().setConstant(
        std::sqrt(2.));

    return ivec;
}
}  // namespace KelvinVector_detail

template <int KelvinVectorSize>
const Eigen::Matrix<double, KelvinVectorSize, KelvinVectorSize>
    Invariants<KelvinVectorSize>::deviatoric_projection =
        KelvinVector_detail::initDeviatoricProjection<KelvinVectorSize>();

template <int KelvinVectorSize>
Eigen::Matrix<double, KelvinVectorSize, KelvinVectorSize> const
    Invariants<KelvinVectorSize>::spherical_projection =
        KelvinVector_detail::initSphericalProjection<KelvinVectorSize>();

template <int KelvinVectorSize>
const Eigen::Matrix<double, KelvinVectorSize, 1>
    Invariants<KelvinVectorSize>::identity2 =
        KelvinVector_detail::initIdentity2<KelvinVectorSize>();

template <int KelvinVectorSize>
const Eigen::Matrix<double, KelvinVectorSize, 1>
    Invariants<KelvinVectorSize>::ones2 =
        KelvinVector_detail::initOnes2<KelvinVectorSize>();

}  // namespace KelvinVector
}  // namespace MathLib
