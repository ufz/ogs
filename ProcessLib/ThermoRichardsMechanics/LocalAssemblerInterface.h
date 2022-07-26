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
            static_cast<int>(this->integration_method_.getIntegrationOrder()))
        {
            OGS_FATAL(
                "Setting integration point initial conditions; The integration "
                "order of the local assembler for element {:d} is different "
                "from the integration order in the initial condition.",
                this->element_.getID());
        }

        if (name == "sigma_ip")
        {
            if (this->process_data_.initial_stress != nullptr)
            {
                OGS_FATAL(
                    "Setting initial conditions for stress from integration "
                    "point data and from a parameter '{:s}' is not possible "
                    "simultaneously.",
                    this->process_data_.initial_stress->name);
            }
            return ProcessLib::setIntegrationPointKelvinVectorData<
                DisplacementDim>(values, this->current_states_,
                                 [](auto const& cs)
                                 { return cs.s_mech_data.sigma_eff; });
        }

        if (name == "saturation_ip")
        {
            return ProcessLib::setIntegrationPointScalarData(
                values, this->current_states_,
                [](auto& state) -> auto& { return state.S_L_data.S_L; });
        }
        if (name == "porosity_ip")
        {
            return ProcessLib::setIntegrationPointScalarData(
                values, this->current_states_,
                [](auto& state) -> auto& { return state.poro_data.phi; });
        }
        if (name == "transport_porosity_ip")
        {
            return ProcessLib::setIntegrationPointScalarData(
                values, this->current_states_, [](auto& state) -> auto& {
                    return state.transport_poro_data.phi;
                });
        }
        if (name == "swelling_stress_ip")
        {
            return ProcessLib::setIntegrationPointKelvinVectorData<
                DisplacementDim>(values, this->current_states_,
                                 [](auto const& cs)
                                 { return cs.swelling_data.sigma_sw; });
        }
        if (name == "epsilon_ip")
        {
            return ProcessLib::setIntegrationPointKelvinVectorData<
                DisplacementDim>(values, this->current_states_,
                                 [](auto const& cs)
                                 { return cs.eps_data.eps; });
        }
        return 0;
    }

    virtual std::vector<double> getSigma() const = 0;

    virtual std::vector<double> const& getIntPtSigma(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> getSwellingStress() const = 0;

    virtual std::vector<double> const& getIntPtSwellingStress(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> getEpsilon() const = 0;

    virtual std::vector<double> const& getIntPtEpsilon(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> getSaturation() const = 0;

    virtual std::vector<double> const& getIntPtSaturation(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> getPorosity() const = 0;

    virtual std::vector<double> const& getIntPtPorosity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> getTransportPorosity() const = 0;

    virtual std::vector<double> const& getIntPtTransportPorosity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDryDensitySolid(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtLiquidDensity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtViscosity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    // TODO move to NumLib::ExtrapolatableElement
    virtual unsigned getNumberOfIntegrationPoints() const = 0;

    virtual typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned /*integration_point*/) const = 0;

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
