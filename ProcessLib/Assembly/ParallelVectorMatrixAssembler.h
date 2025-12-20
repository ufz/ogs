// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "BaseLib/ContainerTools.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/AbstractJacobianAssembler.h"
#include "ProcessLib/Assembly/MatrixOutput.h"

namespace ProcessLib::Assembly
{
class ParallelVectorMatrixAssembler
{
public:
    explicit ParallelVectorMatrixAssembler(
        AbstractJacobianAssembler& jacobian_assembler);

    void assembleWithJacobian(
        BaseLib::PolymorphicRandomAccessContainerView<
            LocalAssemblerInterface> const& local_assemblers,
        std::vector<std::size_t> const& active_elements,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
        const double t, double const dt, std::vector<GlobalVector*> const& xs,
        std::vector<GlobalVector*> const& x_prevs, int const process_id,
        GlobalVector& b, GlobalMatrix& Jac);

private:
    AbstractJacobianAssembler& jacobian_assembler_;
    LocalMatrixOutput local_matrix_output_;
    GlobalMatrixOutput global_matrix_output_;

    int const num_threads_;
};

}  // namespace ProcessLib::Assembly
