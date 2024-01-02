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
        std::vector<
            std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
            dof_tables,
        const double t, double const dt, std::vector<GlobalVector*> const& xs,
        std::vector<GlobalVector*> const& x_prevs, int const process_id,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac);

private:
    AbstractJacobianAssembler& jacobian_assembler_;
    LocalMatrixOutput local_matrix_output_;
    GlobalMatrixOutput global_matrix_output_;

    int const num_threads_;
};

}  // namespace ProcessLib::Assembly
