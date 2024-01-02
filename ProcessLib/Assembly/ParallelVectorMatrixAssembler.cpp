/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ParallelVectorMatrixAssembler.h"

#include <cstdlib>
#include <fstream>

#include "BaseLib/StringTools.h"
#include "BaseLib/ThreadException.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Assembly/MatrixAssemblyStats.h"
#include "ProcessLib/Assembly/MatrixElementCache.h"
#include "ProcessLib/Assembly/MatrixOutput.h"

namespace
{
void assembleWithJacobianOneElement(
    const std::size_t mesh_item_id,
    ProcessLib::LocalAssemblerInterface& local_assembler,
    const NumLib::LocalToGlobalIndexMap& dof_table, const double t,
    const double dt, const GlobalVector& x, const GlobalVector& x_prev,
    std::vector<double>& local_M_data, std::vector<double>& local_K_data,
    std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
    std::vector<GlobalIndexType>& indices,
    ProcessLib::AbstractJacobianAssembler& jacobian_assembler,
    ProcessLib::Assembly::MultiMatrixElementCache& cache)
{
    indices = NumLib::getIndices(mesh_item_id, dof_table);

    local_M_data.clear();
    local_K_data.clear();
    local_b_data.clear();
    local_Jac_data.clear();

    auto const local_x = x.get(indices);
    auto const local_x_prev = x_prev.get(indices);
    jacobian_assembler.assembleWithJacobian(
        local_assembler, t, dt, local_x, local_x_prev, local_M_data,
        local_K_data, local_b_data, local_Jac_data);

    if (local_Jac_data.empty())
    {
        OGS_FATAL(
            "No Jacobian has been assembled! This might be due to "
            "programming errors in the local assembler of the "
            "current process.");
    }

    cache.add(local_M_data, local_K_data, local_b_data, local_Jac_data,
              indices);
}

int getNumberOfThreads()
{
    char const* const num_threads_env = std::getenv("OGS_ASM_THREADS");

    if (!num_threads_env)
    {
        return 1;
    }

    if (std::strlen(num_threads_env) == 0)
    {
        OGS_FATAL("The environment variable OGS_ASM_THREADS is set but empty.");
    }

    std::string num_threads_str{num_threads_env};
    BaseLib::trim(num_threads_str);

    std::istringstream num_threads_iss{num_threads_str};
    int num_threads = -1;

    num_threads_iss >> num_threads;

    if (!num_threads_iss)
    {
        OGS_FATAL("Error parsing OGS_ASM_THREADS (= \"{}\").", num_threads_env);
    }

    if (!num_threads_iss.eof())
    {
        OGS_FATAL(
            "Error parsing OGS_ASM_THREADS (= \"{}\"): not read entirely, the "
            "remainder is \"{}\"",
            num_threads_env,
            num_threads_iss.str().substr(num_threads_iss.tellg()));
    }

    if (num_threads < 1)
    {
        OGS_FATAL(
            "You asked (via OGS_ASM_THREADS) to assemble with {} < 1 thread.",
            num_threads);
    }

    return num_threads;
}
}  // namespace

namespace ProcessLib::Assembly
{
ParallelVectorMatrixAssembler::ParallelVectorMatrixAssembler(
    AbstractJacobianAssembler& jacobian_assembler)
    : jacobian_assembler_{jacobian_assembler},
      num_threads_(getNumberOfThreads())
{
}

void ParallelVectorMatrixAssembler::assembleWithJacobian(
    BaseLib::PolymorphicRandomAccessContainerView<
        LocalAssemblerInterface> const& local_assemblers,
    std::vector<std::size_t> const& active_elements,
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
        dof_tables,
    const double t, double const dt, std::vector<GlobalVector*> const& xs,
    std::vector<GlobalVector*> const& x_prevs, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    // checks //////////////////////////////////////////////////////////////////
    if (process_id != 0)
    {
        OGS_FATAL("Process id is not 0 but {}", process_id);
    }

    if (dof_tables.size() != 1)
    {
        OGS_FATAL("More than 1 dof table");
    }
    auto const& dof_table = dof_tables.front().get();

    if (xs.size() != 1)
    {
        OGS_FATAL("More than 1 solution vector");
    }
    auto const& x = *xs.front();

    if (x_prevs.size() != 1)
    {
        OGS_FATAL("More than 1 x_prev vector");
    }
    auto const& x_prev = *x_prevs.front();

    // algorithm ///////////////////////////////////////////////////////////////

    auto stats = CumulativeStats<MultiStats>::create();
    ConcurrentMatrixView M_view(M);
    ConcurrentMatrixView K_view(K);
    ConcurrentMatrixView b_view(b);
    ConcurrentMatrixView Jac_view(Jac);

    ThreadException exception;
    bool assembly_error = false;
#pragma omp parallel num_threads(num_threads_)
    {
#ifdef _OPENMP
#pragma omp single nowait
        {
            INFO("Number of threads: {}", omp_get_num_threads());
        }
#endif

        // temporary data only stored here in order to avoid frequent memory
        // reallocations.
        std::vector<double> local_M_data;
        std::vector<double> local_K_data;
        std::vector<double> local_b_data;
        std::vector<double> local_Jac_data;
        std::vector<GlobalIndexType> indices;

        // copy to avoid concurrent access
        auto const jac_asm = jacobian_assembler_.copy();
        auto stats_this_thread = stats->clone();

        MultiMatrixElementCache cache{M_view, K_view, b_view, Jac_view,
                                      stats_this_thread->data};

        // TODO corner case: what if all elements on a submesh are deactivated?
        if (active_elements.empty())
        {
            // due to MSVC++ error:
            // error C3016: 'element_id': index variable in OpenMP 'for'
            // statement must have signed integral type
            std::ptrdiff_t const n_loc_asm =
                static_cast<std::ptrdiff_t>(local_assemblers.size());

#pragma omp for nowait
            for (std::ptrdiff_t element_id = 0; element_id < n_loc_asm;
                 ++element_id)
            {
                if (assembly_error)
                {
                    continue;
                }
                auto& loc_asm = local_assemblers[element_id];

                try
                {
                    assembleWithJacobianOneElement(
                        element_id, loc_asm, dof_table, t, dt, x, x_prev,
                        local_M_data, local_K_data, local_b_data,
                        local_Jac_data, indices, *jac_asm, cache);
                }
                catch (...)
                {
                    exception.capture();
                    assembly_error = true;
                    continue;
                }

                local_matrix_output_(t, process_id, element_id, local_M_data,
                                     local_K_data, local_b_data,
                                     &local_Jac_data);
            }
        }
        else
        {
            // due to MSVC++ error:
            // error C3016: 'i': index variable in OpenMP 'for' statement must
            // have signed integral type
            std::ptrdiff_t const n_act_elem =
                static_cast<std::ptrdiff_t>(active_elements.size());

#pragma omp for nowait
            for (std::ptrdiff_t i = 0; i < n_act_elem; ++i)
            {
                if (assembly_error)
                {
                    continue;
                }

                auto const element_id = active_elements[i];
                auto& loc_asm = local_assemblers[element_id];

                try
                {
                    assembleWithJacobianOneElement(
                        element_id, loc_asm, dof_table, t, dt, x, x_prev,
                        local_M_data, local_K_data, local_b_data,
                        local_Jac_data, indices, *jac_asm, cache);
                }
                catch (...)
                {
                    exception.capture();
                    assembly_error = true;
                    continue;
                }

                local_matrix_output_(t, process_id, element_id, local_M_data,
                                     local_K_data, local_b_data,
                                     &local_Jac_data);
            }
        }
    }  // OpenMP parallel section

    stats->print();

    global_matrix_output_(t, process_id, M, K, b, &Jac);
    exception.rethrow();
}
}  // namespace ProcessLib::Assembly
