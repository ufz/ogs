/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/Process.h"

#include "LocalAssemblerInterface.h"
#include "ThermoMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoMechanics
{
struct ThermoMechanicsLocalAssemblerInterface;

struct SigmaIntegrationPointWriter final : public IntegrationPointWriter
{
    explicit SigmaIntegrationPointWriter(
        int const n_components,
        int const integration_order,
        std::function<std::vector<std::vector<double>>()>
            callback)
        : _n_components(n_components),
          _integration_order(integration_order),
          _callback(callback)
    {
    }

    int numberOfComponents() const override { return _n_components; }
    int integrationOrder() const override { return _integration_order; }

    std::string name() const override
    {
        // TODO (naumov) remove ip suffix. Probably needs modification of the
        // mesh properties, s.t. there is no "overlapping" with cell/point data.
        // See getOrCreateMeshProperty.
        return "sigma_ip";
    }

    std::vector<std::vector<double>> values() const override
    {
        return _callback();
    }

private:
    int const _n_components;
    int const _integration_order;
    std::function<std::vector<std::vector<double>>()> _callback;
};

template <int DisplacementDim>
class ThermoMechanicsProcess final : public Process
{
public:
    ThermoMechanicsProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        ThermoMechanicsProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller,
        bool const use_monolithic_scheme);

    //! \name ODESystem interface
    //! @{
    bool isLinear() const override;
    //! @}

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(
        const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
        GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

    void preTimestepConcreteProcess(
        GlobalVector const& x, double const t, double const dt,
        const int process_id) override;

    void postTimestepConcreteProcess(GlobalVector const& x, const double t,
                                     const double delta_t,
                                     int const process_id) override;

private:
    std::vector<MeshLib::Node*> _base_nodes;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_base_nodes;
    ThermoMechanicsProcessData<DisplacementDim> _process_data;

    std::vector<std::unique_ptr<ThermoMechanicsLocalAssemblerInterface>>
        _local_assemblers;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;
};

extern template class ThermoMechanicsProcess<2>;
extern template class ThermoMechanicsProcess<3>;

}  // namespace ThermoMechanics
}  // namespace ProcessLib
