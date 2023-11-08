/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 6, 2023, 1:42 PM
 */

#pragma once

#include "MechanicsBase.h"
#include "ParameterLib/Parameter.h"

namespace MaterialLib
{
namespace Solids
{
/**
 * \brief Linear transverse isotropic elastic model.
 *
 *   The parameters of the linear transverse isotropic elastic model are
 *    - \f$E_{i}\f$, the Young’s modulus within the plane of isotropy,
 *    - \f$E_{a}\f$, the Young’s modulus w.r.t. the direction of anisotropy,
 *    - \f$\nu_{ii}\f$, the Poisson’s ratio within the plane of isotropy,
 *    - \f$\nu_{ia}\f$, the Poisson ratio perpendicular to the plane of
 * isotropy, due to strain in the plane of isotropy,
 *    - \f$G_{ia}\f$, the shear modulus between directions of isotropy and
 * anisotropy.
 *
 *    where the subscript \f$i\f$ means isotropy, and the subscript \f$a\f$
 * means anisotropy.
 *
 *   With the given parameter, the in-plane shear modulus, \f$G_{ii}\f$ is
 * computed as
 *    \f[
 *      G_{ii} = \frac{E_{i}}{2(1+\nu_{ii})},
 *    \f]
 *    while the  in-plane Poisson ratio, \f$\nu_{ai}\f$, which is due to the
 * strain perpendicular to the plane of isotropy, is calculated by the following
 * equation:
 *     \f[
 *        \nu_{ai} = \nu_{ia} \frac{E_{a}}{E_{i}}.
 *     \f]
 *
 *  For 3D problems, assuming the plane of isotropy to be spanned by the basis
 *  vectors \f$\mathbf{e}_1\f$ and \f$\mathbf{e}_2\f$, respectively, and the
 *  direction of anisotropy is defined by the basis vector \f$\mathbf{e}_3\f$,
 *  the following relations hold:
 *    \f{eqnarray*}{
 *        E_{i} &=& E_1 & =& E_2, \\
 *        E_{a} &=& E_3, & &\\
 *        \nu_{ii} &=& \nu_{12} & =& \nu_{21}, \\
 *        \nu_{ia} &=& \nu_{13} & =& \nu_{23}, \\
 *        \nu_{ai} &=& \nu_{31} & =& \nu_{32}, \\
 *        G_{ia} &=& G_{13} & =& G_{23},\\
 *        G_{ai} &=& G_{ia}. & &
 *    \f}
 *    Under such assumption, the matrix form of the
 *    elastic tensor for strain and stress in the Kelvin vector in the local
 *    system is
 *    \f[
 *      \begin{bmatrix}
 *          a_{ii} & b_{ii} &b_{ai} & 0 & 0 & 0\\
 *          b_{ii} & a_{ii} &b_{ai} & 0 & 0 & 0\\
 *          b_{ai} & b_{ai} &a_{ai} & 0 & 0 & 0\\
 *          0 & 0  & 0 & 2 c_{ii} & 0 & 0\\
 *          0 & 0  & 0 & 0 & 2 c_{ai} & 0\\
 *          0 & 0  & 0 & 0 & 0 & 2 c_{ai}
 *      \end{bmatrix}.
 *    \f]
 *    The matrix elements are:
 *    \f{eqnarray*}{
 *        a_{ii} &=& \frac{1-\nu_{ia}\nu_{ai}}{E_{i} E_{a} D}, \\
 *        a_{ai} &=& \frac{1-\nu_{ii}^2}{E_{i}^2 D}, \\
 *        b_{ii} &=& \frac{\nu_{ii}+\nu_{ia}\nu_{ai}}{E_{i} E_{a} D}, \\
 *        b_{ai} &=& \frac{\nu_{ia}(1+\nu_{ii})}{E_{i}^2 D}, \\
 *        c_{ii} &=& \frac{E_{i}}{2(1+\nu_{ii})}, \\
 *        c_{ai} &=& G_{ia},
 *    \f}
 *   with
 *   \f[
 *        D = \frac{(1+\nu_{ii})(1-\nu_{ii}-2\nu_{ia}\nu_{ai})}{E_{i}^2E_{a}}.
 *   \f]
 *   (also see Chapter 9.1 in \cite kolditz2012thermo).
 *
 *   For plane strain problems, assuming the direction of anisotropy is defined
 *   by the basis vector \f$\mathbf{e}_1\f$, the plane of isotropy to be spanned
 *   by the basis vector \f$\mathbf{e}_0\f$ and the unit off-plane direction
 *   \f$\mathbf{e}_2\f$, the
 *    following relations hold:
 *    \f{eqnarray*}{
 *        E_{i} &=& E_1 & =& E_3, \\
 *        E_{a} &=& E_2, & &\\
 *        \nu_{ii} &=& \nu_{13} & =& \nu_{31}, \\
 *        \nu_{ia} &=& \nu_{32} & =& \nu_{12}, \\
 *        \nu_{ai} &=& \nu_{23} & =& \nu_{21}, \\
 *        G_{ia} &=& G_{32} & =& G_{12},\\
 *        G_{ai} &=& G_{ia}. & &
 *    \f}
 *
 *   Based on this assumption, the matrix form of the elastic tensor for strain
 *   and stress in the Kelvin vector in the local system is
 *    \f[
 *      \begin{bmatrix}
 *          a_{ii} & b_{ai} &b_{ii} & 0\\
 *          b_{ai} & a_{ai} &b_{ai} & 0\\
 *          b_{ii} & b_{ai} &a_{ii} & 0\\
 *          0 & 0  & 0 & 2 c_{ai}
 *      \end{bmatrix}
 *    \f]
 *  for plane strain problems.
 *
 *  \note This model requires the definition of a local coordinate system. A
 * local coordinate system can be defined with explicit or implicit bases. The
 * unit norm to the transverse isotropy plane must be input as the parameter
 * value for the base \"basis_vector_1\" for a 2D local coordinate system, or
 * for the base \"basis_vector_2\" for a 3D local coordinate system.
 *  @see ParameterLib#CoordinateSystem.
 */
template <int DisplacementDim>
class LinearElasticTransverseIsotropic : public MechanicsBase<DisplacementDim>
{
public:
    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    using P = ParameterLib::Parameter<double>;

    LinearElasticTransverseIsotropic(
        P const& E_i, P const& E_a, P const& nu_ii, P const& nu_ia,
        P const& G_ia,
        std::optional<ParameterLib::CoordinateSystem> const&
            local_coordinate_system)
        : E_i_p_(E_i),
          E_a_p_(E_a),
          nu_ii_p_(nu_ii),
          nu_ia_p_(nu_ia),
          G_ia_p_(G_ia),
          local_coordinate_system_(local_coordinate_system)
    {
    }

    double computeFreeEnergyDensity(
        double const /*t*/,
        ParameterLib::SpatialPosition const& /*x*/,
        double const /*dt*/,
        KelvinVector const& eps,
        KelvinVector const& sigma,
        typename MechanicsBase<DisplacementDim>::
            MaterialStateVariables const& /* material_state_variables */)
        const override
    {
        return eps.dot(sigma) / 2;
    }

    std::optional<
        std::tuple<typename MechanicsBase<DisplacementDim>::KelvinVector,
                   std::unique_ptr<typename MechanicsBase<
                       DisplacementDim>::MaterialStateVariables>,
                   typename MechanicsBase<DisplacementDim>::KelvinMatrix>>
    integrateStress(
        MaterialPropertyLib::VariableArray const& variable_array_prev,
        MaterialPropertyLib::VariableArray const& variable_array,
        double const t, ParameterLib::SpatialPosition const& x, double const dt,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const override;

    KelvinMatrix getElasticTensor(double const t,
                                  ParameterLib::SpatialPosition const& x,
                                  double const T) const;

    double getBulkModulus(double const t,
                          ParameterLib::SpatialPosition const& x,
                          KelvinMatrix const* const /*C*/) const override
    {
        const double E_i = E_i_p_(t, x)[0];
        const double E_a = E_a_p_(t, x)[0];
        const double nu_i = nu_ii_p_(t, x)[0];
        const double nu_ia = nu_ia_p_(t, x)[0];

        // Average Young's modulus
        double const E_av = 2. * E_i / 3. + E_a / 3.;

        // Poisson ratio in the plane of isotropy, due to the strain
        // perpendicular to the plane of isotropy.  nu_ai=nu_ia*E_a/E_i
        const double nu_ai = nu_ia * E_a / E_i;

        // Average Poisson ratio
        //    12     13    21   23   31    32
        //    ai     ai    ia   ii   ia    ii
        const double nu_av = (nu_ai + nu_ia + nu_i) / 3.0;

        return E_av / 3 / (1 - 2 * nu_av);
    }

protected:
    /// It is the in-plane Young’s modulus, \f$E_{i}\f$.
    P const& E_i_p_;
    /// It is the Young’s modulus w.r.t. the direction of anisotropy,
    /// \f$E_{a}\f$.
    P const& E_a_p_;
    /// It is the in-plane Poisson’s ratio, \f$\nu_{ii}\f$.
    P const& nu_ii_p_;
    /// It is the Poisson ratio perpendicular to the plane of isotropy, due to
    /// strain in the plane of isotropy, \f$\nu_{ia}\f$.
    P const& nu_ia_p_;
    /// It is the shear modulus between directions of isotropy and anisotropy,
    /// \f$G_{ia}\f$.
    P const& G_ia_p_;
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system_;

private:
    KelvinMatrix getElasticTensorLeftTopCorner(
        double const t, ParameterLib::SpatialPosition const& x) const;
};

}  // namespace Solids
}  // namespace MaterialLib
