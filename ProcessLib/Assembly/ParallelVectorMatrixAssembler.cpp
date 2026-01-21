// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ParallelVectorMatrixAssembler.h"

#include <cstdlib>
#include <fstream>
#include <range/v3/view/iota.hpp>
#include <vector>

#include "BaseLib/OgsAsmThreads.h"
#include "BaseLib/StringTools.h"
#include "BaseLib/ThreadException.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Assembly/MatrixAssemblyStats.h"
#include "ProcessLib/Assembly/MatrixElementCache.h"
#include "ProcessLib/Assembly/MatrixOutput.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"
#include "ProcessLib/LocalAssemblerInterface.h"

namespace
{
void assembleOneElement(const std::size_t mesh_item_id,
                        ProcessLib::LocalAssemblerInterface& local_assembler,
                        const NumLib::LocalToGlobalIndexMap& dof_table,
                        const double t, const double dt, const GlobalVector& x,
                        const GlobalVector& x_prev,
                        std::vector<double>& local_M_data,
                        std::vector<double>& local_K_data,
                        std::vector<double>& local_b_data,
                        ProcessLib::Assembly::MultiMatrixElementCache<2>& cache)
{
    std::vector<GlobalIndexType> indices =
        NumLib::getIndices(mesh_item_id, dof_table);

    local_M_data.clear();
    local_K_data.clear();
    local_b_data.clear();

    auto const local_x = x.get(indices);
    auto const local_x_prev = x_prev.get(indices);
    local_assembler.assemble(t, dt, local_x, local_x_prev, local_M_data,
                             local_K_data, local_b_data);

    cache.add(local_M_data, local_K_data, local_b_data, std::move(indices));
}

void assembleForStaggeredSchemeOneElement(
    const std::size_t mesh_item_id,
    ProcessLib::LocalAssemblerInterface& local_assembler,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
    const double t, const double dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    std::vector<double>& local_M_data, std::vector<double>& local_K_data,
    std::vector<double>& local_b_data,
    ProcessLib::Assembly::MultiMatrixElementCache<2>& cache)
{
    std::vector<std::vector<GlobalIndexType>> indices_of_processes;
    indices_of_processes.reserve(dof_tables.size());
    transform(cbegin(dof_tables), cend(dof_tables),
              back_inserter(indices_of_processes), [&](auto const* dof_table)
              { return NumLib::getIndices(mesh_item_id, *dof_table); });

    auto local_coupled_xs =
        ProcessLib::getCoupledLocalSolutions(x, indices_of_processes);
    auto const local_x = MathLib::toVector(local_coupled_xs);

    auto local_coupled_x_prevs =
        ProcessLib::getCoupledLocalSolutions(x_prev, indices_of_processes);
    auto const local_x_prev = MathLib::toVector(local_coupled_x_prevs);

    local_M_data.clear();
    local_K_data.clear();
    local_b_data.clear();

    local_assembler.assembleForStaggeredScheme(t, dt, local_x, local_x_prev,
                                               process_id, local_M_data,
                                               local_K_data, local_b_data);

    auto const& indices = indices_of_processes[process_id];
    cache.add(local_M_data, local_K_data, local_b_data, indices);
}

void assembleWithJacobianOneElement(
    const std::size_t mesh_item_id,
    ProcessLib::LocalAssemblerInterface& local_assembler,
    const NumLib::LocalToGlobalIndexMap& dof_table, const double t,
    const double dt, const GlobalVector& x, const GlobalVector& x_prev,
    std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
    ProcessLib::AbstractJacobianAssembler& jacobian_assembler,
    ProcessLib::Assembly::MultiMatrixElementCache<1>& cache)
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

void assembleWithJacobianForStaggeredSchemeOneElement(
    const std::size_t mesh_item_id,
    ProcessLib::LocalAssemblerInterface& local_assembler,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
    const double t, const double dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
    ProcessLib::AbstractJacobianAssembler& jacobian_assembler,
    ProcessLib::Assembly::MultiMatrixElementCache<1>& cache)
{
    std::vector<std::vector<GlobalIndexType>> indices_of_processes;
    indices_of_processes.reserve(dof_tables.size());
    transform(cbegin(dof_tables), cend(dof_tables),
              back_inserter(indices_of_processes),
              [&](auto const* dof_table)
              { return NumLib::getIndices(mesh_item_id, *dof_table); });

    auto local_coupled_xs =
        ProcessLib::getCoupledLocalSolutions(x, indices_of_processes);
    auto const local_x = MathLib::toVector(local_coupled_xs);

    auto local_coupled_x_prevs =
        ProcessLib::getCoupledLocalSolutions(x_prev, indices_of_processes);
    auto const local_x_prev = MathLib::toVector(local_coupled_x_prevs);

    std::vector<GlobalIndexType> const& indices =
        indices_of_processes[process_id];

    local_b_data.clear();
    local_Jac_data.clear();

    jacobian_assembler.assembleWithJacobianForStaggeredScheme(
        local_assembler, t, dt, local_x, local_x_prev, process_id, local_b_data,
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

//! Runs the passed \c assemble functor for each active local assembler.
void runAssemblyImpl(
    BaseLib::PolymorphicRandomAccessContainerView<
        ProcessLib::LocalAssemblerInterface> const& local_assemblers,
    ThreadException& exception,
    auto local_matrix_output,
    auto assemble,
    std::ptrdiff_t const num_active_local_assemblers,
    auto get_element_id_callback)
{
#pragma omp for nowait
    for (std::ptrdiff_t i = 0; i < num_active_local_assemblers; ++i)
    {
        [[unlikely]] if (exception)
        {
            continue;
        }

        auto const element_id = get_element_id_callback(i);
        auto& loc_asm = local_assemblers[element_id];

        try
        {
            assemble(element_id, loc_asm);
        }
        catch (...)
        {
            exception.capture();
            continue;
        }

        local_matrix_output(element_id);
    }
}

//! Runs the passed \c assemble functor for each active local assembler.
//!
//! Forwards all calls to runAssemblyImpl().
void runAssembly(
    BaseLib::PolymorphicRandomAccessContainerView<
        ProcessLib::LocalAssemblerInterface> const& local_assemblers,
    std::vector<std::size_t> const* const active_elements,
    ThreadException& exception,
    auto local_matrix_output,
    auto assemble)
{
    [[likely]] if (!active_elements)
    {
        std::ptrdiff_t const n_elements =
            static_cast<std::ptrdiff_t>(local_assemblers.size());

        runAssemblyImpl(local_assemblers, exception, local_matrix_output,
                        assemble, n_elements, [](auto i) { return i; });
    }
    else
    {
        std::ptrdiff_t const n_elements =
            static_cast<std::ptrdiff_t>(active_elements->size());

        runAssemblyImpl(local_assemblers, exception, local_matrix_output,
                        assemble, n_elements, [active_elements](auto i)
                        { return (*active_elements)[i]; });
    }
}
}  // namespace

namespace ProcessLib::Assembly
{
ParallelVectorMatrixAssembler::ParallelVectorMatrixAssembler(
    AbstractJacobianAssembler& jacobian_assembler)
    : jacobian_assembler_{jacobian_assembler},
      num_threads_(BaseLib::getNumberOfAssemblyThreads())
{
    INFO("Threads used for ParallelVectorMatrixAssembler: {}.", num_threads_);
}

void ParallelVectorMatrixAssembler::assemble(
    BaseLib::PolymorphicRandomAccessContainerView<
        LocalAssemblerInterface> const& local_assemblers,
    std::vector<std::size_t> const* const active_elements,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
    const double t, double const dt, std::vector<GlobalVector*> const& xs,
    std::vector<GlobalVector*> const& x_prevs, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    // checks //////////////////////////////////////////////////////////////////
    if (dof_tables.size() != xs.size())
    {
        OGS_FATAL("Different number of DOF tables and solution vectors.");
    }

    std::size_t const number_of_processes = xs.size();

    // algorithm ///////////////////////////////////////////////////////////////

    auto stats = CumulativeStats<MultiStats<2>>::create(num_threads_);

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
        std::vector<double> local_M_data;
        std::vector<double> local_K_data;
        std::vector<double> local_b_data;

        // copy to avoid concurrent access
        auto const jac_asm = jacobian_assembler_.copy();

        auto stats_this_thread = stats->clone();
        MultiMatrixElementCache<2> cache{M, K, b, stats_this_thread->data,
                                         num_threads_};

        auto local_matrix_output = [&](std::ptrdiff_t element_id)
        {
            local_matrix_output_(t, process_id, element_id, local_M_data,
                                 local_K_data, local_b_data);
        };

        // Monolithic scheme
        if (number_of_processes == 1)
        {
            assert(process_id == 0);
            auto const& dof_table = *dof_tables[0];
            auto const& x = *xs[0];
            auto const& x_prev = *x_prevs[0];

            runAssembly(local_assemblers, active_elements, exception,
                        local_matrix_output,
                        [&](auto element_id, auto& loc_asm)
                        {
                            assembleOneElement(element_id, loc_asm, dof_table,
                                               t, dt, x, x_prev, local_M_data,
                                               local_K_data, local_b_data,
                                               cache);
                        });
        }
        else  // Staggered scheme
        {
            runAssembly(local_assemblers, active_elements, exception,
                        local_matrix_output,
                        [&](auto element_id, auto& loc_asm)
                        {
                            assembleForStaggeredSchemeOneElement(
                                element_id, loc_asm, dof_tables, t, dt, xs,
                                x_prevs, process_id, local_M_data, local_K_data,
                                local_b_data, cache);
                        });
        }
    }

    global_matrix_output_(t, process_id, M, K, b);
    exception.rethrow();
}

void ParallelVectorMatrixAssembler::assembleWithJacobian(
    BaseLib::PolymorphicRandomAccessContainerView<
        LocalAssemblerInterface> const& local_assemblers,
    std::vector<std::size_t> const* const active_elements,
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

    std::size_t const number_of_processes = xs.size();
    // algorithm ///////////////////////////////////////////////////////////////

    auto stats = CumulativeStats<MultiStats<1>>::create(num_threads_);

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

        MultiMatrixElementCache<1> cache{b, Jac, stats_this_thread->data,
                                         num_threads_};

        auto local_matrix_output = [&](std::ptrdiff_t element_id)
        {
            local_matrix_output_(t, process_id, element_id, local_b_data,
                                 local_Jac_data);
        };

        // Monolithic scheme
        if (number_of_processes == 1)
        {
            assert(process_id == 0);
            auto const& dof_table = *dof_tables[0];
            auto const& x = *xs[0];
            auto const& x_prev = *x_prevs[0];

            runAssembly(
                local_assemblers, active_elements, exception,
                local_matrix_output,
                [&](auto element_id, auto& loc_asm)
                {
                    assembleWithJacobianOneElement(
                        element_id, loc_asm, dof_table, t, dt, x, x_prev,
                        local_b_data, local_Jac_data, *jac_asm, cache);
                });
        }
        else  // Staggered scheme
        {
            runAssembly(
                local_assemblers, active_elements, exception,
                local_matrix_output,
                [&](auto element_id, auto& loc_asm)
                {
                    assembleWithJacobianForStaggeredSchemeOneElement(
                        element_id, loc_asm, dof_tables, t, dt, xs, x_prevs,
                        process_id, local_b_data, local_Jac_data, *jac_asm,
                        cache);
                });
        }
    }

    stats->print();

    global_matrix_output_(t, process_id, b, Jac);
    exception.rethrow();
}
}  // namespace ProcessLib::Assembly
