/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "LocalAssemblerInterface.h"
#include "ProcessLib/Process.h"
#include "SmallDeformationProcessData.h"

namespace ProcessLib
{
namespace SmallDeformation
{
template <int DisplacementDim>
class SmallDeformationProcess final : public Process
{
public:
    SmallDeformationProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        SmallDeformationProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables);

    //! \name ODESystem interface
    //! @{
    bool isLinear() const override;
    //! @}

private:
    using LocalAssemblerInterface =
        SmallDeformationLocalAssemblerInterface<DisplacementDim>;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(double const t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& x_prev,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalVector& b, GlobalMatrix& Jac) override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     std::vector<GlobalVector*> const& x_prev,
                                     double const t, double const dt,
                                     int const process_id) override;

    void computeSecondaryVariableConcrete(double const t, double const dt,
                                          std::vector<GlobalVector*> const& x,
                                          GlobalVector const& x_prev,
                                          int const process_id) override;

private:
    SmallDeformationProcessData<DisplacementDim> _process_data;

    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    MeshLib::PropertyVector<double>* _nodal_forces = nullptr;
    MeshLib::PropertyVector<double>* _material_forces = nullptr;

    Assembly::GlobalMatrixOutput _global_output;
};

extern template class SmallDeformationProcess<2>;
extern template class SmallDeformationProcess<3>;

}  // namespace SmallDeformation
}  // namespace ProcessLib
