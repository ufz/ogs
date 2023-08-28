/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "AssembledMatrixCache.h"

#include "BaseLib/RunTime.h"
#include "MathLib/LinAlg/FinalizeMatrixAssembly.h"
#include "MathLib/LinAlg/FinalizeVectorAssembly.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"

namespace ProcessLib::ComponentTransport
{

AssembledMatrixCache::AssembledMatrixCache(bool const is_linear,
                                           bool const use_monolithic_scheme)
    : is_linear_{is_linear && use_monolithic_scheme}
{
    if (is_linear && !use_monolithic_scheme)
    {
        OGS_FATAL(
            "You requested to assemble only once in combination with staggered "
            "coupling. This use case is not yet implemented.");
    }

    if (is_linear_)
    {
        WARN(
            "You specified that the ComponentTransport process is linear. With "
            "that optimization the process will be assembled only once and the "
            "non-linear solver will do only one iteration per time step. No "
            "non-linearities will be resolved and OGS will not detect if there "
            "are any non-linearities. It is your responsibility to ensure that "
            "the assembled equation systems are linear, indeed! There is no "
            "safety net!");
    }
}

void AssembledMatrixCache::assemble(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
        dof_tables,
    VectorMatrixAssembler& global_assembler,
    std::vector<
        std::unique_ptr<ComponentTransportLocalAssemblerInterface>> const&
        local_assemblers,
    std::vector<std::size_t> const& active_element_ids)
{
    if (bool const cache_empty = K_ == nullptr; cache_empty)
    {
        BaseLib::RunTime time_asm;
        time_asm.start();

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeSelectedMemberDereferenced(
            global_assembler, &VectorMatrixAssembler::assemble,
            local_assemblers, active_element_ids, dof_tables, t, dt, x, x_prev,
            process_id, M, K, b);

        MathLib::finalizeMatrixAssembly(M);
        MathLib::finalizeMatrixAssembly(K);
        MathLib::finalizeVectorAssembly(b);

        INFO("[time] Calling local assemblers took {:g} s", time_asm.elapsed());

        if (is_linear_)
        {
            DBUG("Saving global K, M, b for later reuse.");

            BaseLib::RunTime time_save;
            time_save.start();

            K_ = MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(K);
            M_ = MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(M);
            b_ = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(b);

            INFO("[time] Saving global K, M, b took {:g} s",
                 time_save.elapsed());
        }
    }
    else
    {
        DBUG("Reusing saved global K, M, b.");

        BaseLib::RunTime time_restore;
        time_restore.start();

        MathLib::LinAlg::copy(*K_, K);
        MathLib::LinAlg::copy(*M_, M);
        MathLib::LinAlg::copy(*b_, b);

        INFO("[time] Restoring global K, M, b took {:g} s",
             time_restore.elapsed());
    }
}
}  // namespace ProcessLib::ComponentTransport
