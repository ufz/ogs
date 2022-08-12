/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveSetting.h"
#include "ProcessLib/ThermoRichardsMechanics/ThermoRichardsMechanicsProcessData.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "ProcessLib/Utils/TransposeInPlace.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct LocalAssemblerInterface : public ProcessLib::LocalAssemblerInterface,
                                 public NumLib::ExtrapolatableElement
{
    LocalAssemblerInterface(
        MeshLib::Element const& e,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        ThermoRichardsMechanicsProcessData<DisplacementDim>& process_data)
        : process_data_(process_data),
          integration_method_(integration_method),
          element_(e),
          is_axially_symmetric_(is_axially_symmetric),
          solid_material_(MaterialLib::Solids::selectSolidConstitutiveRelation(
              process_data_.solid_materials, process_data_.material_ids,
              e.getID()))
    {
        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();

        current_states_.resize(n_integration_points);
        prev_states_.resize(n_integration_points);
        output_data_.resize(n_integration_points);

        material_states_.reserve(n_integration_points);
        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            material_states_.emplace_back(solid_material_);
        }
    }

    std::size_t setIPDataInitialConditions(std::string const& name,
                                           double const* values,
                                           int const integration_order)
    {
        if (integration_order !=
            static_cast<int>(integration_method_.getIntegrationOrder()))
        {
            OGS_FATAL(
                "Setting integration point initial conditions; The integration "
                "order of the local assembler for element {:d} is different "
                "from the integration order in the initial condition.",
                element_.getID());
        }

        if (name == "sigma_ip")
        {
            if (process_data_.initial_stress != nullptr)
            {
                OGS_FATAL(
                    "Setting initial conditions for stress from integration "
                    "point data and from a parameter '{:s}' is not possible "
                    "simultaneously.",
                    process_data_.initial_stress->name);
            }
            return ProcessLib::setIntegrationPointKelvinVectorData<
                DisplacementDim>(
                values, current_states_,
                [](auto& cs) -> auto& { return cs.s_mech_data.sigma_eff; });
        }

        if (name == "saturation_ip")
        {
            return ProcessLib::setIntegrationPointScalarData(
                values, current_states_,
                [](auto& state) -> auto& { return state.S_L_data.S_L; });
        }
        if (name == "porosity_ip")
        {
            return ProcessLib::setIntegrationPointScalarData(
                values, current_states_,
                [](auto& state) -> auto& { return state.poro_data.phi; });
        }
        if (name == "transport_porosity_ip")
        {
            return ProcessLib::setIntegrationPointScalarData(
                values, current_states_, [](auto& state) -> auto& {
                    return state.transport_poro_data.phi;
                });
        }
        if (name == "swelling_stress_ip")
        {
            return ProcessLib::setIntegrationPointKelvinVectorData<
                DisplacementDim>(
                values, current_states_,
                [](auto& cs) -> auto& { return cs.swelling_data.sigma_sw; });
        }
        if (name == "epsilon_ip")
        {
            return ProcessLib::setIntegrationPointKelvinVectorData<
                DisplacementDim>(
                values, current_states_,
                [](auto& cs) -> auto& { return cs.eps_data.eps; });
        }
        return 0;
    }

private:
    std::vector<double> getSigma() const
    {
        constexpr int kelvin_vector_size =
            MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

        return transposeInPlace<kelvin_vector_size>(
            [this](std::vector<double>& values)
            { return getIntPtSigma(0, {}, {}, values); });
    }

    std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            current_states_, [](auto const& cs) -> auto const& {
                return cs.s_mech_data.sigma_eff;
            },
            cache);
    }

    std::vector<double> getSwellingStress() const
    {
        constexpr int kelvin_vector_size =
            MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

        return transposeInPlace<kelvin_vector_size>(
            [this](std::vector<double>& values)
            { return getIntPtSwellingStress(0, {}, {}, values); });
    }

    std::vector<double> const& getIntPtSwellingStress(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
    {
        constexpr int kelvin_vector_size =
            MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
        auto const n_integration_points = current_states_.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, n_integration_points);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            auto const& sigma_sw = current_states_[ip].swelling_data.sigma_sw;
            cache_mat.col(ip) =
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_sw);
        }

        return cache;
    }

    std::vector<double> getEpsilon() const
    {
        constexpr int kelvin_vector_size =
            MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

        return transposeInPlace<kelvin_vector_size>(
            [this](std::vector<double>& values)
            { return getIntPtEpsilon(0, {}, {}, values); });
    }

    std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            current_states_,
            [](auto const& cs) -> auto const& { return cs.eps_data.eps; },
            cache);
    }

    std::vector<double> const& getIntPtDarcyVelocity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
    {
        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();

        cache.clear();
        auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, DisplacementDim, n_integration_points);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            cache_matrix.col(ip).noalias() =
                output_data_[ip].darcy_data.v_darcy;
        }

        return cache;
    }

    std::vector<double> getSaturation() const
    {
        std::vector<double> result;
        getIntPtSaturation(0, {}, {}, result);
        return result;
    }

    std::vector<double> const& getIntPtSaturation(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
    {
        return ProcessLib::getIntegrationPointScalarData(
            current_states_,
            [](auto const& state) -> auto const& { return state.S_L_data.S_L; },
            cache);
    }

    std::vector<double> getPorosity() const
    {
        std::vector<double> result;
        getIntPtPorosity(0, {}, {}, result);
        return result;
    }

    std::vector<double> const& getIntPtPorosity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
    {
        return ProcessLib::getIntegrationPointScalarData(
            current_states_, [](auto const& state) -> auto const& {
                return state.poro_data.phi;
            },
            cache);
    }

    std::vector<double> getTransportPorosity() const
    {
        std::vector<double> result;
        getIntPtTransportPorosity(0, {}, {}, result);
        return result;
    }

    std::vector<double> const& getIntPtTransportPorosity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
    {
        return ProcessLib::getIntegrationPointScalarData(
            current_states_,
            [](auto const& state) -> auto const& {
                return state.transport_poro_data.phi;
            },
            cache);
    }

    std::vector<double> const& getIntPtDryDensitySolid(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
    {
        return ProcessLib::getIntegrationPointScalarData(
            output_data_,
            [](auto const& out) -> auto const& {
                return out.rho_S_data.dry_density_solid;
            },
            cache);
    }

    std::vector<double> const& getIntPtLiquidDensity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
    {
        return ProcessLib::getIntegrationPointScalarData(
            output_data_,
            [](auto const& out) -> auto const& {
                return out.rho_L_data.rho_LR;
            },
            cache);
    }

    std::vector<double> const& getIntPtViscosity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
    {
        return ProcessLib::getIntegrationPointScalarData(
            output_data_, [](auto const& out) -> auto const& {
                return out.mu_L_data.viscosity;
            },
            cache);
    }

public:
    // TODO move to NumLib::ExtrapolatableElement
    unsigned getNumberOfIntegrationPoints() const
    {
        return integration_method_.getNumberOfPoints();
    }

    typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned integration_point) const
    {
        return *material_states_[integration_point].material_state_variables;
    }

    void postTimestepConcrete(Eigen::VectorXd const& /*local_x*/,
                              double const /*t*/,
                              double const /*dt*/) override
    {
        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            // TODO re-evaluate part of the assembly in order to be consistent?
            material_states_[ip].pushBackState();
        }

        prev_states_ = current_states_;
    }

    struct IPDataAccessorForExtrapolation
    {
        std::string name;
        unsigned num_comp;
        std::function<std::vector<double> const&(
            LocalAssemblerInterface<DisplacementDim> const&,
            const double /*t*/,
            std::vector<GlobalVector*> const& /*x*/,
            std::vector<
                NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
            std::vector<double>& /*cache*/)>
            ip_data_accessor;
    };

    static std::vector<IPDataAccessorForExtrapolation>
    getIPDataAccessorsForExtrapolation()
    {
        using Self = LocalAssemblerInterface<DisplacementDim>;
        constexpr auto kv_size =
            MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

        return {{"sigma", kv_size, &Self::getIntPtSigma},
                {"swelling_stress", kv_size, &Self::getIntPtSwellingStress},
                {"epsilon", kv_size, &Self::getIntPtEpsilon},
                {"velocity", DisplacementDim, &Self::getIntPtDarcyVelocity},
                {"saturation", 1, &Self::getIntPtSaturation},
                {"porosity", 1, &Self::getIntPtPorosity},
                {"transport_porosity", 1, &Self::getIntPtTransportPorosity},
                {"dry_density_solid", 1, &Self::getIntPtDryDensitySolid},
                {"liquid_density", 1, &Self::getIntPtLiquidDensity},
                {"viscosity", 1, &Self::getIntPtViscosity}};
    }

    struct IPDataAccessorForIPWriter
    {
        std::string name;
        unsigned num_comp;
        std::vector<double> (LocalAssemblerInterface<DisplacementDim>::*
                                 ip_data_accessor)() const;
    };

    static std::vector<IPDataAccessorForIPWriter>
    getIPDataAccessorsForIPWriter()
    {
        using Self = LocalAssemblerInterface<DisplacementDim>;
        constexpr auto kv_size =
            MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

        // TODO (naumov) remove ip suffix. Probably needs modification of the
        // mesh properties, s.t. there is no "overlapping" with cell/point data.
        // See getOrCreateMeshProperty.
        return {
            {"sigma_ip", kv_size, &Self::getSigma},
            {"saturation_ip", 1, &Self::getSaturation},
            {"porosity_ip", 1, &Self::getPorosity},
            {"transport_porosity_ip", 1, &Self::getTransportPorosity},
            {"swelling_stress_ip", kv_size, &Self::getSwellingStress},
            {"epsilon_ip", kv_size, &Self::getEpsilon},
        };
    }

protected:
    ThermoRichardsMechanicsProcessData<DisplacementDim>& process_data_;

    std::vector<StatefulData<DisplacementDim>>
        current_states_;  // TODO maybe do not store but rather re-evaluate for
                          // state update
    std::vector<StatefulData<DisplacementDim>> prev_states_;

    // Material state is special, because it contains both the current and the
    // old state.
    std::vector<MaterialStateData<DisplacementDim>> material_states_;

    NumLib::GenericIntegrationMethod const& integration_method_;
    MeshLib::Element const& element_;
    bool const is_axially_symmetric_;

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material_;

    std::vector<OutputData<DisplacementDim>> output_data_;
};

}  // namespace ProcessLib::ThermoRichardsMechanics
