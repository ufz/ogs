/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "LargeDeformationProcessData.h"
#include "LocalAssemblerInterface.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace LargeDeformation
{
template <int DisplacementDim>
class LargeDeformationProcess final : public Process
{
public:
    LargeDeformationProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        LargeDeformationProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables);

    //! \name ODESystem interface
    //! @{
    bool isLinear() const override;
    //! @}

private:
    using LocalAssemblerInterface =
        LargeDeformationLocalAssemblerInterface<DisplacementDim>;

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
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac) override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     std::vector<GlobalVector*> const& x_prev,
                                     double const t, double const dt,
                                     int const process_id) override;

    void computeSecondaryVariableConcrete(double const t, double const dt,
                                          std::vector<GlobalVector*> const& x,
                                          GlobalVector const& x_prev,
                                          int const process_id) override;

private:
    LargeDeformationProcessData<DisplacementDim> _process_data;

    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    MeshLib::PropertyVector<double>* _nodal_forces = nullptr;
    MeshLib::PropertyVector<double>* _material_forces = nullptr;

    Assembly::GlobalMatrixOutput _global_output;
};

extern template class LargeDeformationProcess<2>;
extern template class LargeDeformationProcess<3>;

}  // namespace LargeDeformation
}  // namespace ProcessLib
