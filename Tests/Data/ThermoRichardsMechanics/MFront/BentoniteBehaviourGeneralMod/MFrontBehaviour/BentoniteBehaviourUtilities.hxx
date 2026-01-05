// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

/**
 * \file
 * This header defines a set of auxiliary functions for the conversion of the
 * gradients, the thermodynamic forces and the consistent tangent operator.
 */

#ifndef LIB_TFEL_MATERIAL_BENTONITEBEHAVIOURUTILITIES_HXX
#define LIB_TFEL_MATERIAL_BENTONITEBEHAVIOURUTILITIES_HXX

#include "TFEL/Math/stensor.hxx"
#include "TFEL/Math/tvector.hxx"

namespace bentonite_behaviour_utilities
{

//! a simple alias
using real = double;
//! a simple alias
using size_type = unsigned short;

namespace internals
{

template <size_type N>
struct NormalizeTangentOperator;

template <>
struct NormalizeTangentOperator<2u>
{
    template <typename Stensor4>
    static void exe(Stensor4& Dt)
    {
        constexpr const auto cste = tfel::math::Cste<real>::sqrt2;
        Dt(0, 3) *= cste;
        Dt(1, 3) *= cste;
        Dt(2, 3) *= cste;
        Dt(3, 0) *= cste;
        Dt(3, 1) *= cste;
        Dt(3, 2) *= cste;
        Dt(3, 3) *= 2;
    }  // end of exe
};  // end of struct NormalizeTangentOperator<2u>

template <>
struct NormalizeTangentOperator<3u>
{
    template <typename Stensor4>
    static void exe(Stensor4& Dt)
    {
        constexpr const auto cste = tfel::math::Cste<real>::sqrt2;
        Dt(0, 3) *= cste;
        Dt(1, 3) *= cste;
        Dt(2, 3) *= cste;
        Dt(0, 4) *= cste;
        Dt(1, 4) *= cste;
        Dt(2, 4) *= cste;
        Dt(0, 5) *= cste;
        Dt(1, 5) *= cste;
        Dt(2, 5) *= cste;
        Dt(3, 0) *= cste;
        Dt(3, 1) *= cste;
        Dt(3, 2) *= cste;
        Dt(4, 0) *= cste;
        Dt(4, 1) *= cste;
        Dt(4, 2) *= cste;
        Dt(5, 0) *= cste;
        Dt(5, 1) *= cste;
        Dt(5, 2) *= cste;
        Dt(3, 3) *= 2;
        Dt(3, 4) *= 2;
        Dt(3, 5) *= 2;
        Dt(4, 3) *= 2;
        Dt(4, 4) *= 2;
        Dt(4, 5) *= 2;
        Dt(5, 3) *= 2;
        Dt(5, 4) *= 2;
        Dt(5, 5) *= 2;
    }  // end of exe
};  // end of struct NormalizeTangentOperator<2u>

}  // end of namespace internals

/*!
 * \brief convert the gradient from `MFront` and store them in an array
 * that is to be passed to the Bentonite behaviour.
 * This array will contain:
 * [eps11, eps22, eps33, gamma12, gamma13, gamma23, s, T]
 * \note One must take into account that `MFront` stores the strain as follows:
 * [eps11, eps22, eps33, sqrt(2) * eps12, sqrt(2) * eps13, sqrt(2) * eps23]
 * \note The temperature is not a gradient in the `MFront` sense
 * \param[in] e: strain
 * \param[in] p_L: liquid pressure
 * \param[in] T: temperature
 * \return an array containing the gradients for the Bentonite behaviour
 */
template <size_type N>
tfel::math::tvector<8u, real> convert_gradients_from_mfront(
    const tfel::math::stensor<N, real>& e, const real p_L, const real T)
{
    constexpr auto size = tfel::math::StensorDimeToSize<N>::value;
    constexpr const auto cste = tfel::math::Cste<real>::sqrt2;
    auto r = tfel::math::tvector<8u, real>{};
    tfel::fsalgo::copy<size>::exe(e.begin(), r.begin());
    tfel::fsalgo::fill<6u - size>::exe(r.begin() + e.size(), real(0));
    for (size_type i = 3; i != e.size(); ++i)
    {
        r[i] *= cste;
    }
    r[6] = p_L / 1000.0; /* converting the the ogs pressure in Pa into kPa
                            required by the MFront/generalmod */
    r[7] = T;
    return r;
}  // end of convert_gradients_from_mfront

/*!
 * \brief convert the thermodynamic forces from `MFront` and store them in an
 * array that is to be passed to the Bentonite behaviour. This array will
 * contain: [signet11, signet22, signet33, signet12, signet13, signet23, Sr]
 * \note One must take into account that `MFront` stores the stress as follows:
 * [sig11, sig22, sig33, sqrt(2) * sig12, sqrt(2) * sig13, sqrt(2) * sig23]
 * \param[in] sig: stress
 * \param[in] Sr: saturation
 * \return an array containing the thermodynamic forces for the Bentonite
 * behaviour
 */
template <size_type N>
tfel::math::tvector<7u, real> convert_thermodynamic_forces_from_mfront(
    const tfel::math::stensor<N, real>& sig, const real Sr)
{
    constexpr auto size = tfel::math::StensorDimeToSize<N>::value;
    constexpr const auto icste = tfel::math::Cste<real>::isqrt2;
    auto r = tfel::math::tvector<7u, real>{};
    tfel::fsalgo::copy<size>::exe(sig.begin(), r.begin());
    tfel::fsalgo::fill<6u - size>::exe(r.begin() + sig.size(), real(0));
    for (size_type i = 3; i != sig.size(); ++i)
    {
        r[i] *= icste;
    }
    for (size_type i = 0; i != sig.size(); ++i)
    {
        r[i] /= 1000.0; /* Pa to kPa conversion */
    }
    r[6] = Sr;
    return r;
}  // end of convert_gradients_from_mfront

/*!
 * \brief convert the thermodynamic forces computed for the Bentonite behaviour
 * to `MFront`. \param[out] sig: stress \param[out] Sr: saturation \param[in]
 * tf: the thermodynamic forces computed for the Bentonite behaviour
 */
template <size_type N>
void convert_thermodynamic_forces_to_mfront(
    tfel::math::stensor<N, real>& sig,
    real& Sr,
    const tfel::math::tvector<7u, real>& tf)
{
    constexpr auto size = tfel::math::StensorDimeToSize<N>::value;
    constexpr const auto cste = tfel::math::Cste<real>::sqrt2;
    tfel::fsalgo::copy<size>::exe(tf.begin(), sig.begin());
    for (size_type i = 3; i != sig.size(); ++i)
    {
        sig[i] *= cste;
    }
    for (size_type i = 0; i != sig.size(); ++i)
    {
        sig[i] *= 1000.0; /* kPa to Pa conversion */
    }
    Sr = tf[6];
}  // end of convert_gradients_from_mfront

/*!
 * \brief convert the derivative of the stress tensor with respect to the strain
 * tensor to the `MFront` conventions. \tparam N: space dimension \tparam
 * Stensor4: type of the fourth order tensor \param[in,out] Dt: fourth order
 * tensor to be converted
 */
template <size_type N, typename Stensor4>
void normalize_tangent_operator(Stensor4& Dt)
{
    internals::NormalizeTangentOperator<N>::exe(Dt);
}  // end of normalize_tangent_operator

/*!
 * \brief convert the tangent operator returned by the Bentonite behaviour to
 * the tangent operator blocks defined in `MFront` \tparam Stensor4: type of the
 * fourth order tensor \tparam Stensor: type of the second order tensor
 * \param[out] dsig_ddeto: derivative of the stress tensor with respect to the
 * strain tensor \param[out] dsig_dds: derivative of the stress tensor with
 * respect to the suction \param[out] dSr_dds: derivative of the saturation with
 * respect to the suction \param[in] Dtg: tangent operator returned by the
 * Bentonite behaviour
 */
template <typename Stensor4, typename Stensor>
void convert_tangent_operator(Stensor4& dsig_ddeto,
                              Stensor& dsig_ddp_L,
                              real& dSr_ddp_L,
                              const tfel::math::tmatrix<7, 8, real>& Dtg)
{
    constexpr const auto cste = tfel::math::Cste<real>::sqrt2;
    constexpr const auto N = tfel::math::StensorTraits<Stensor>::dime;
    for (size_type i = 0; i != tfel::math::StensorDimeToSize<N>::value; ++i)
    {
        for (size_type j = 0; j != tfel::math::StensorDimeToSize<N>::value; ++j)
        {
            dsig_ddeto(i, j) = 1000.0 * Dtg(i, j); /* kPa to Pa coversion */
        }
        dsig_ddp_L[i] = Dtg(i, 6);
    }
    normalize_tangent_operator<N>(dsig_ddeto);
    // normalise the derivative of the stress tensor with respect to the suction
    for (size_type i = 3; i != dsig_ddp_L.size(); ++i)
    {
        dsig_ddp_L[i] *= cste;
    }
    dSr_ddp_L =
        Dtg(6, 6) /
        1000.0; /* suction to liquid pressure, and kPa to Pa conversion */
}  // end of convert_tangent_operator

}  // end of namespace bentonite_behaviour_utilities

#endif /* LIB_TFEL_MATERIAL_BENTONITEBEHAVIOURUTILITIES_HXX */
