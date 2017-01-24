/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <logog/include/logog.hpp>

#include "BaseLib/Error.h"

#include "KelvinVector.h"
#include "MechanicsBase.h"

namespace MaterialLib
{
namespace Solids
{
//
// Variables specific to the material model.
//
struct Lubby2MaterialProperties
{
    using P = ProcessLib::Parameter<double>;
    Lubby2MaterialProperties(P const& GK0_,
                             P const& GM0_,
                             P const& KM0_,
                             P const& etaK0_,
                             P const& etaM0_,
                             P const& mK_,
                             P const& mvK_,
                             P const& mvM_)
        : GK0(GK0_),
          GM0(GM0_),
          KM0(KM0_),
          etaK0(etaK0_),
          etaM0(etaM0_),
          mK(mK_),
          mvK(mvK_),
          mvM(mvM_)
    {
    }

    // basic material parameters
    P const& GK0;
    P const& GM0;
    P const& KM0;
    P const& etaK0;
    P const& etaM0;
    P const& mK;
    P const& mvK;
    P const& mvM;
};

namespace detail
{
template <int DisplacementDim>
struct LocalLubby2Properties
{
    LocalLubby2Properties(double const t,
                          ProcessLib::SpatialPosition const& x,
                          Lubby2MaterialProperties const& mp)
        : GM0(mp.GM0(t, x)[0]),
          KM0(mp.KM0(t, x)[0]),
          GK0(mp.GK0(t, x)[0]),
          etaK0(mp.etaK0(t, x)[0]),
          etaM0(mp.etaM0(t, x)[0]),
          mK(mp.mK(t, x)[0]),
          mvK(mp.mvK(t, x)[0]),
          mvM(mp.mvM(t, x)[0])
    {
    }

    void update(double const s_eff)
    {
        GK = GK0 * std::exp(mK * GM0 * s_eff);
        etaK = etaK0 * std::exp(mvK * GM0 * s_eff);
        etaM = etaM0 * std::exp(mvM * GM0 * s_eff);
    }

    double const GM0;
    double const KM0;
    double const GK0;
    double const etaK0;
    double const etaM0;
    double const mK;
    double const mvK;
    double const mvM;

    // Solution dependent values.
    double GK;
    double etaK;
    double etaM;
};
}

template <int DisplacementDim>
class Lubby2 final : public MechanicsBase<DisplacementDim>
{
public:
    struct MaterialStateVariables
        : public MechanicsBase<DisplacementDim>::MaterialStateVariables
    {
        MaterialStateVariables()
        {
            eps_K_t.resize(KelvinVectorSize);
            eps_K_j.resize(KelvinVectorSize);
            eps_M_t.resize(KelvinVectorSize);
            eps_M_j.resize(KelvinVectorSize);
        }

        void loadState()
        {
            eps_K_j = eps_K_t;
            eps_M_j = eps_M_t;
        }

        void pushBackState() override
        {
            eps_K_t = eps_K_j;
            eps_M_t = eps_M_j;
        }

        using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;
        using KelvinMatrix = ProcessLib::KelvinMatrixType<DisplacementDim>;
        /// Deviatoric strain in the viscous kelvin element during the current
        /// iteration
        KelvinVector eps_K_t;
        KelvinVector eps_K_j;
        /// Deviatoric strain in the viscous maxwell element during the current
        /// iteration
        KelvinVector eps_M_t;
        KelvinVector eps_M_j;
    };

    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() override
    {
        return std::unique_ptr<
            typename MechanicsBase<DisplacementDim>::MaterialStateVariables>{
            new MaterialStateVariables};
    }

public:
    static int const KelvinVectorSize =
        ProcessLib::KelvinVectorDimensions<DisplacementDim>::value;
    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix = ProcessLib::KelvinMatrixType<DisplacementDim>;

    static int const JacobianResidualSize =
        3 * KelvinVectorSize;  // Three is the number of components in the
                               // jacobian/residual, not the space dimension.
    using ResidualVector = Eigen::Matrix<double, JacobianResidualSize, 1>;
    using JacobianMatrix = Eigen::Matrix<double,
                                         JacobianResidualSize,
                                         JacobianResidualSize,
                                         Eigen::RowMajor>;

public:
    explicit Lubby2(Lubby2MaterialProperties& material_properties)
        : _mp(material_properties)
    {
    }

    bool computeConstitutiveRelation(
        double const t,
        ProcessLib::SpatialPosition const& x_position,
        double const dt,
        KelvinVector const& eps_prev,
        KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        KelvinVector& sigma,
        KelvinMatrix& C,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
            material_state_variables) override;

private:
    /// Calculates the 18x1 residual vector.
    void calculateResidualBurgers(
        double const dt,
        const KelvinVector& strain_curr,
        const KelvinVector& stress_curr,
        KelvinVector& strain_Kel_curr,
        const KelvinVector& strain_Kel_t,
        KelvinVector& strain_Max_curr,
        const KelvinVector& strain_Max_t,
        ResidualVector& res,
        detail::LocalLubby2Properties<DisplacementDim> const& properties);

    /// Calculates the 18x18 Jacobian.
    void calculateJacobianBurgers(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        JacobianMatrix& Jac,
        double s_eff,
        const KelvinVector& sig_i,
        const KelvinVector& eps_K_i,
        detail::LocalLubby2Properties<DisplacementDim> const& properties);

    /// Calculates the 18x6 derivative of the residuals with respect to total
    /// strain.
    ///
    /// Function definition can not be moved into implementation because of a
    /// MSVC compiler errors. See
    /// http://stackoverflow.com/questions/1484885/strange-vc-compile-error-c2244
    /// and https://support.microsoft.com/en-us/kb/930198
    Eigen::
        Matrix<double, JacobianResidualSize, KelvinVectorSize, Eigen::RowMajor>
        calculatedGdEBurgers() const
    {
        Eigen::Matrix<double,
                      JacobianResidualSize,
                      KelvinVectorSize,
                      Eigen::RowMajor>
            dGdE;
        dGdE.setZero();
        dGdE.template block<KelvinVectorSize, KelvinVectorSize>(0, 0)
            .diagonal()
            .setConstant(-2);
        return dGdE;
    }

private:
    Lubby2MaterialProperties _mp;
};

extern template class Lubby2<2>;
extern template class Lubby2<3>;

}  // namespace Solids
}  // namespace MaterialLib

#include "Lubby2-impl.h"
