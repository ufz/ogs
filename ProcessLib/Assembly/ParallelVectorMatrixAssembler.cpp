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
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>
#include <vector>

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
    std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
    ProcessLib::AbstractJacobianAssembler& jacobian_assembler,
    ProcessLib::Assembly::MultiMatrixElementCache& cache)
{
    std::vector<GlobalIndexType> const& indices =
        NumLib::getIndices(mesh_item_id, dof_table);

    local_b_data.clear();
    local_Jac_data.clear();

    auto const local_x = x.get(indices);
    auto const local_x_prev = x_prev.get(indices);
    jacobian_assembler.assembleWithJacobian(local_assembler, t, dt, local_x,
                                            local_x_prev, local_b_data,
                                            local_Jac_data);

    if (local_Jac_data.empty())
    {
        OGS_FATAL(
            "No Jacobian has been assembled! This might be due to "
            "programming errors in the local assembler of the "
            "current process.");
    }

    cache.add(local_b_data, local_Jac_data, indices);
}

/// Returns a vector of active element ids with corresponding local assemblers.
std::vector<
    std::tuple<std::ptrdiff_t,
               std::reference_wrapper<ProcessLib::LocalAssemblerInterface>>>
collectActiveLocalAssemblers(
    BaseLib::PolymorphicRandomAccessContainerView<
        ProcessLib::LocalAssemblerInterface> const& local_assemblers,
    std::vector<std::size_t> const& active_elements)
{
    auto id_and_local_asm = [&local_assemblers](std::size_t const id)
        -> std::tuple<std::ptrdiff_t, std::reference_wrapper<
                                          ProcessLib::LocalAssemblerInterface>>
    { return {id, local_assemblers[id]}; };

    auto create_ids_asm_pairs = [&](auto const& element_ids)
    {
        return element_ids | ranges::views::transform(id_and_local_asm) |
               ranges::to<std::vector>();
    };

    if (active_elements.empty())
    {
        return create_ids_asm_pairs(ranges::views::iota(
            static_cast<std::size_t>(0), local_assemblers.size()));
    }
    return create_ids_asm_pairs(active_elements);
}

void runAssemblyForEachLocalAssembler(
    std::vector<std::tuple<
        std::ptrdiff_t,
        std::reference_wrapper<ProcessLib::LocalAssemblerInterface>>> const&
        ids_local_assemblers,
    NumLib::LocalToGlobalIndexMap const& dof_table, double const t,
    double const dt, GlobalVector const& x, GlobalVector const& x_prev,
    std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
    ProcessLib::AbstractJacobianAssembler& jac_asm, ThreadException& exception,
    ProcessLib::Assembly::MultiMatrixElementCache& cache,
    auto local_matrix_output)
{
    // due to MSVC++ error:
    // error C3016: 'element_id': index variable in OpenMP 'for'
    // statement must have signed integral type
    std::ptrdiff_t n_elements =
        static_cast<std::ptrdiff_t>(ids_local_assemblers.size());
#pragma omp for nowait
    for (std::ptrdiff_t i = 0; i < n_elements; ++i)
    {
        if (exception)
        {
            continue;
        }
        auto [element_id, loc_asm] = ids_local_assemblers[i];

        try
        {
            assembleWithJacobianOneElement(element_id, loc_asm, dof_table, t,
                                           dt, x, x_prev, local_b_data,
                                           local_Jac_data, jac_asm, cache);
        }
        catch (...)
        {
            exception.capture();
            continue;
        }

        local_matrix_output(element_id);
    }
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
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
    const double t, double const dt, std::vector<GlobalVector*> const& xs,
    std::vector<GlobalVector*> const& x_prevs, int const process_id,
    GlobalVector& b, GlobalMatrix& Jac)
{
    // checks //////////////////////////////////////////////////////////////////
    if (dof_tables.size() != xs.size())
    {
        OGS_FATAL("Different number of DOF tables and solution vectors.");
    }

    auto const& dof_table = *dof_tables[process_id];
    auto const& x = *xs[process_id];
    auto const& x_prev = *x_prevs[process_id];

    // algorithm ///////////////////////////////////////////////////////////////

    auto stats = CumulativeStats<MultiStats>::create();
    ConcurrentMatrixView b_view(b);
    ConcurrentMatrixView Jac_view(Jac);

    ThreadException exception;
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
        std::vector<double> local_b_data;
        std::vector<double> local_Jac_data;

        // copy to avoid concurrent access
        auto const jac_asm = jacobian_assembler_.copy();
        auto stats_this_thread = stats->clone();

        MultiMatrixElementCache cache{b_view, Jac_view,
                                      stats_this_thread->data};

        auto local_matrix_output = [&](std::ptrdiff_t element_id)
        {
            local_matrix_output_(t, process_id, element_id, local_b_data,
                                 local_Jac_data);
        };

        // TODO corner case: what if all elements on a submesh are deactivated?
        runAssemblyForEachLocalAssembler(
            collectActiveLocalAssemblers(local_assemblers, active_elements),
            dof_table, t, dt, x, x_prev, local_b_data, local_Jac_data, *jac_asm,
            exception, cache, local_matrix_output);
    }

    stats->print();

    global_matrix_output_(t, process_id, b, Jac);
    exception.rethrow();
}
}  // namespace ProcessLib::Assembly
