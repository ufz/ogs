/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <limits>

#include "MechanicsBase.h"
#include "ParameterLib/Parameter.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
class LinearElasticOrthotropic : public MechanicsBase<DisplacementDim>
{
public:
    struct EvaluatedMaterialProperties
    {
        /// Youngs moduli for 1-based indices.
        constexpr double E(int const i) const
        {
            assert(i == 1 || i == 2 || i == 3);
            return youngs_moduli[i - 1];
        }

        /// Shear moduli for 1-based indices.
        constexpr double G(int const i, int const j) const
        {
            assert(i == 1 || i == 2 || i == 3);
            assert(j == 1 || j == 2 || j == 3);
            if (i == 1 && j == 2)
            {
                return shear_moduli[0];
            }
            if (i == 2 && j == 3)
            {
                return shear_moduli[1];
            }
            if (i == 1 && j == 3)
            {
                return shear_moduli[2];
            }

            return std::numeric_limits<double>::quiet_NaN();
        }

        /// Poisson's ratios for 1-based indices.
        constexpr double nu(int const i, int const j) const
        {
            assert(i == 1 || i == 2 || i == 3);
            assert(j == 1 || j == 2 || j == 3);
            if (i == 1 && j == 2)
            {
                return poissons_ratios[0];
            }
            if (i == 2 && j == 3)
            {
                return poissons_ratios[1];
            }
            if (i == 1 && j == 3)
            {
                return poissons_ratios[2];
            }

            // The next set is scaled by Youngs modulus s.t.
            // nu_ji = nu_ij * E_j/E_i.
            if (i == 2 && j == 1)
            {
                return poissons_ratios[0] * E(i) / E(j);
            }
            if (i == 3 && j == 2)
            {
                return poissons_ratios[1] * E(i) / E(j);
            }
            if (i == 3 && j == 1)
            {
                return poissons_ratios[2] * E(i) / E(j);
            }

            return std::numeric_limits<double>::quiet_NaN();
        }

        std::array<double, 3> youngs_moduli;
        std::array<double, 3> shear_moduli;
        std::array<double, 3> poissons_ratios;
    };

    /// Variables specific to the material model
    struct MaterialProperties
    {
        using P = ParameterLib::Parameter<double>;
        using X = ParameterLib::SpatialPosition;

        MaterialProperties(P const& youngs_moduli_, P const& shear_moduli_,
                           P const& poissons_ratios_)
            : youngs_moduli(youngs_moduli_),
              shear_moduli(shear_moduli_),
              poissons_ratios(poissons_ratios_)
        {
        }

        EvaluatedMaterialProperties evaluate(
            double const t, ParameterLib::SpatialPosition const& x) const
        {
            auto const E = youngs_moduli(t, x);
            auto const G = shear_moduli(t, x);
            auto const nu = poissons_ratios(t, x);
            return {
                {E[0], E[1], E[2]}, {G[0], G[1], G[2]}, {nu[0], nu[1], nu[2]}};
        }

        P const& youngs_moduli;
        P const& shear_moduli;
        P const& poissons_ratios;  // Stored as nu_12, nu_23, nu_13
    };

public:
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    LinearElasticOrthotropic(
        MaterialProperties material_properties,
        boost::optional<ParameterLib::CoordinateSystem> const&
            local_coordinate_system)
        : _mp(std::move(material_properties)),
          _local_coordinate_system(local_coordinate_system)
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
        double const t, ParameterLib::SpatialPosition const& x,
        double const /*dt*/,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const override;

    KelvinMatrix getElasticTensor(double const t,
                                  ParameterLib::SpatialPosition const& x,
                                  double const T) const;

    MaterialProperties getMaterialProperties() const { return _mp; }

    double getBulkModulus(double const t,
                          ParameterLib::SpatialPosition const& x,
                          KelvinMatrix const* const /*C*/) const override
    {
        auto const& mp = _mp.evaluate(t, x);
        auto const E = [&mp](int const i) { return mp.E(i); };
        auto const nu = [&mp](int const i, int const j) { return mp.nu(i, j); };
        // corresponds to 1/(I:S:I) --> Reuss bound.
        // Voigt bound would proceed as (I:C:I)/9. Use MFront model if you
        // prefer that.
        return E(1) * E(2) * E(3) /
               (E(1) * E(2) + E(1) * E(3) * (1 - 2 * nu(2, 3)) +
                E(2) * E(3) * (1 - 2 * nu(1, 2) * nu(1, 3)));
    }

protected:
    MaterialProperties _mp;
    boost::optional<ParameterLib::CoordinateSystem> const&
        _local_coordinate_system;
};

extern template class LinearElasticOrthotropic<2>;
extern template class LinearElasticOrthotropic<3>;

}  // namespace Solids
}  // namespace MaterialLib
