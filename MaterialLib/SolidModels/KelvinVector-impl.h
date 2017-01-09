/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

namespace MaterialLib
{
namespace SolidModels
{
template <int KelvinVectorSize>
double Invariants<KelvinVectorSize>::equivalentStress(
    Eigen::Matrix<double, KelvinVectorSize, 1> const& deviatoric_v)
{
    assert(std::abs(trace(deviatoric_v)) <
           std::numeric_limits<double>::epsilon());
    return std::sqrt(3 * J2(deviatoric_v));
}

template <int KelvinVectorSize>
double Invariants<KelvinVectorSize>::J2(
    Eigen::Matrix<double, KelvinVectorSize, 1> const& deviatoric_v)
{
    assert(std::abs(trace(deviatoric_v)) <
           std::numeric_limits<double>::epsilon());
    return 0.5 * deviatoric_v.transpose() * deviatoric_v;
}

/// Third invariant, equal to determinant of a deviatoric tensor.
/// \note The input vector must have trace equal zero.
template <int KelvinVectorSize>
double Invariants<KelvinVectorSize>::J3(
    Eigen::Matrix<double, KelvinVectorSize, 1> const& deviatoric_v)
{
    assert(std::abs(trace(deviatoric_v)) <
           std::numeric_limits<double>::epsilon());
    return determinant(deviatoric_v);
}

/// Trace of the corresponding tensor.
template <int KelvinVectorSize>
double Invariants<KelvinVectorSize>::trace(
    Eigen::Matrix<double, KelvinVectorSize, 1> const& v)
{
    return v.template topLeftCorner<3, 1>().sum();
}

//
// Initialization of static Invariant variables.
//

namespace detail
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
}  // namespace detail

template <int KelvinVectorSize>
const Eigen::Matrix<double, KelvinVectorSize, KelvinVectorSize>
    Invariants<KelvinVectorSize>::deviatoric_projection =
        detail::initDeviatoricProjection<KelvinVectorSize>();

template <int KelvinVectorSize>
Eigen::Matrix<double, KelvinVectorSize, KelvinVectorSize> const
    Invariants<KelvinVectorSize>::spherical_projection =
        detail::initSphericalProjection<KelvinVectorSize>();

template <int KelvinVectorSize>
const Eigen::Matrix<double, KelvinVectorSize, 1> Invariants<
    KelvinVectorSize>::identity2 = detail::initIdentity2<KelvinVectorSize>();

}  // namespace SolidModels
}  // namespace MaterialLib
