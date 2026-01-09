// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "BaseLib/MPI.h"
#include "BaseLib/RunTime.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/Exceptions.h"
#include "NumLib/ODESolver/Types.h"
#include "ProcessLib/Assembly/AssembledMatrixCache.h"
#include "ProcessLib/Assembly/AssemblyData.h"
#include "ProcessLib/Assembly/ParallelVectorMatrixAssembler.h"
#include "ProcessLib/ProcessVariable.h"
#include "ProcessLib/Utils/ComputeResiduum.h"
#include "VectorMatrixAssembler.h"

namespace ProcessLib
{
/// Provides basic functionality for the AssemblyMixin that does not require
/// template parameters.
class AssemblyMixinBase
{
    enum class ActiveElementIDsState
    {
        UNINITIALIZED,
        HAS_DEACTIVATED_SUBDOMAINS,
        NO_DEACTIVATED_SUBDOMAINS
    };

protected:
    explicit AssemblyMixinBase(AbstractJacobianAssembler& jacobian_assembler,
                               bool const is_linear,
                               bool const use_monolithic_scheme);

    void initializeAssemblyOnSubmeshes(
        MeshLib::Mesh& bulk_mesh,
        std::vector<std::reference_wrapper<MeshLib::Mesh>> const& submeshes,
        std::vector<std::vector<std::string>> const& residuum_names,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>> const&
            pvs);

    void updateActiveElements(ProcessLib::Process const& process);

    bool isLinear() const { return is_linear_; }

    static void copyResiduumVectorsToBulkMesh(
        GlobalVector const& rhs,
        NumLib::LocalToGlobalIndexMap const& local_to_global_index_map,
        std::vector<std::reference_wrapper<MeshLib::PropertyVector<double>>>
            residuum_vectors);

    static void copyResiduumVectorsToSubmesh(
        int const process_id,
        GlobalVector const& rhs,
        NumLib::LocalToGlobalIndexMap const& local_to_global_index_map,
        ProcessLib::Assembly::SubmeshAssemblyData const& sad);

private:
    void updateActiveElementsImpl(Process const& process);

protected:
    /// SubmeshAssemblyData for each submesh.
    std::vector<ProcessLib::Assembly::SubmeshAssemblyData>
        submesh_assembly_data_;

    /// AssembledMatrixCache for each submesh.
    std::vector<AssembledMatrixCache> submesh_matrix_cache_;

    /// Empty if submesh assembly is used.
    std::optional<ProcessLib::Assembly::BulkMeshAssemblyData>
        bulk_mesh_assembly_data_;

    AssembledMatrixCache bulk_mesh_matrix_cache_;

    Assembly::ParallelVectorMatrixAssembler pvma_;

    /// If the process/problem being simulated is linear, the assembled global
    /// matrices and vectors can be cached and reused.
    bool is_linear_;

    std::optional<NumLib::NonlinearSolverTag> last_assembly_was_;

private:
    ActiveElementIDsState ids_state_ = ActiveElementIDsState::UNINITIALIZED;
};

/**
 * A mixin providing assembly functionality to a specific \c Process.
 *
 * The process must be derived from this class (CRTP).
 */
template <typename Process>
class AssemblyMixin : private AssemblyMixinBase
{
    // Enforce correct use of CRTP, i.e., that Process is derived from
    // AssemblyMixin<Process>.
    using AssemblyMixinBase::AssemblyMixinBase;
    friend Process;

public:
    /**
     * Specifies that the assembly of the process should take place on the
     * passed \c submeshes.
     *
     * \note This method must be called at most once.
     *
     * \attention The passed submeshes must be non-overlapping and must cover
     * the entire simulation domain. The caller is responsible to meet this
     * requirement.
     *
     * \attention The order of the passed \c residuum_names matters and must
     * match the order of the residuum vectors in the process. Otherwise
     * residuum output on submeshes will be wrong.
     *
     * \note As of January 2023 assembly on submeshes has neither been tested
     * with domain decomposition, nor with domain deactivation nor staggered
     * coupling. Use at your own risk in these cases.
     *
     * \todo Implement and test with domain decomposition, domain deactivation
     * and staggered coupling.
     */
    void initializeAssemblyOnSubmeshes(
        std::vector<std::reference_wrapper<MeshLib::Mesh>> const& submeshes,
        std::vector<std::vector<std::string>> const& residuum_names)
    {
        AssemblyMixinBase::initializeAssemblyOnSubmeshes(
            derived().getMesh(), submeshes, residuum_names,
            derived().getProcessVariables());
    }

    void updateActiveElements()
    {
        AssemblyMixinBase::updateActiveElements(derived());
    }

    // sorted_element_subset can be used to restrict assembly to a custom
    // element subset. That's used at the moment (Jan 2026) for some assembly
    // optimizations in the HeatTransportBHE process.
    void assemble(
        double const t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        std::vector<std::size_t> const* const sorted_element_subset = nullptr,
        bool const copy_residua_to_mesh = false)
    {
        DBUG("AssemblyMixin assemble(t={:g}, dt={:g}, process_id={}).", t, dt,
             process_id);

        last_assembly_was_ = NumLib::NonlinearSolverTag::Picard;

        std::vector<NumLib::LocalToGlobalIndexMap const*> const dof_tables =
            derived().getDOFTables(x.size());

        std::exception_ptr const exception =
            submesh_assembly_data_.empty()
                ? assembleOnBulkMesh(t, dt, x, x_prev, dof_tables, process_id,
                                     M, K, b, sorted_element_subset)
                : assembleOnSubmeshes(t, dt, x, x_prev, dof_tables, process_id,
                                      M, K, b, sorted_element_subset,
                                      copy_residua_to_mesh);

        [[unlikely]] if (BaseLib::MPI::anyOf(exception != nullptr))
        {
            if (exception)  // Only the rank with the exception rethrows...
            {
                std::rethrow_exception(exception);
            }
            // TODO should other ranks rather throw a new exception?
            // ... but all ranks quit.
            return;
        }

        if (copy_residua_to_mesh && bulk_mesh_assembly_data_)
        {
            // Note: computeResiduum() uses bwd/fwd Euler scheme.
            GlobalVector const neg_res_bulkmesh = ProcessLib::computeResiduum(
                dt, *x[process_id], *x_prev[process_id], M, K, b);

            AssemblyMixinBase::copyResiduumVectorsToBulkMesh(
                neg_res_bulkmesh, *(dof_tables[process_id]),
                bulk_mesh_assembly_data_->residuum_vectors[process_id]);
        }
    }

    // sorted_element_subset can be used to restrict assembly to a custom
    // element subset. That's used at the moment (Jan 2026) for some assembly
    // optimizations in the HeatTransportBHE process.
    void assembleWithJacobian(
        double const t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalVector& b, GlobalMatrix& Jac,
        std::vector<std::size_t> const* const sorted_element_subset = nullptr,
        bool const copy_residua_to_mesh = false)
    {
        DBUG(
            "AssemblyMixin assembleWithJacobian(t={:g}, dt={:g}, "
            "process_id={}).",
            t, dt, process_id);

        last_assembly_was_ = NumLib::NonlinearSolverTag::Newton;

        // Note: We do not cache the Jacobian or b vector. It does not make
        // sense, because even in the linear case the b vector changes if
        // the solution changes (unless steady state is reached).
        // Todo: If we separated Jacobian and residuum assembly, we could
        // cache at least the Jacobian and only update the residuum.
        if (is_linear_)
        {
            DBUG(
                "You specified that this process is linear. The assembly for "
                "the Newton-Raphson method will be run anyways.");
        }

        std::vector<NumLib::LocalToGlobalIndexMap const*> const dof_tables =
            derived().getDOFTables(x.size());

        std::exception_ptr const exception =
            submesh_assembly_data_.empty()
                ? assembleWithJacobianOnBulkMeshOrOnSubmeshCommon(
                      *bulk_mesh_assembly_data_, t, dt, x, x_prev, dof_tables,
                      process_id, b, Jac, sorted_element_subset)
                : assembleWithJacobianOnSubmeshes(
                      t, dt, x, x_prev, dof_tables, process_id, b, Jac,
                      sorted_element_subset, copy_residua_to_mesh);

        [[unlikely]] if (BaseLib::MPI::anyOf(exception != nullptr))
        {
            if (exception)  // Only the rank with the exception rethrows...
            {
                std::rethrow_exception(exception);
            }
            // ... but all ranks quit.
            return;
        }

        MathLib::LinAlg::finalizeAssembly(Jac);

        // TODO better for submesh case, too?
        // TODO check copy_residua_to_mesh, too. But that needs a refactoring of
        // postTimestep(), computeSecondaryVariable() etc. hooks.
        if (/*copy_residua_to_mesh &&*/ bulk_mesh_assembly_data_)
        {
            AssemblyMixinBase::copyResiduumVectorsToBulkMesh(
                b, *(dof_tables[process_id]),
                bulk_mesh_assembly_data_->residuum_vectors[process_id]);
        }
    }

    void preOutput(double const t,
                   double const dt,
                   std::vector<GlobalVector*> const& x,
                   std::vector<GlobalVector*> const& x_prev,
                   int const process_id)
    {
        if (!last_assembly_was_)
        {
            WARN("AssemblyMixin. Skipping preOutput() step.");
            return;
        }

        switch (*last_assembly_was_)
        {
            case NumLib::NonlinearSolverTag::Picard:
                preOutputPicard(t, dt, x, x_prev, process_id);
                return;
            case NumLib::NonlinearSolverTag::Newton:
                preOutputNewton(t, dt, x, x_prev, process_id);
                return;
        }
    }

private:
    Process& derived() { return static_cast<Process&>(*this); }
    Process const& derived() const
    {
        return static_cast<Process const&>(*this);
    }

    //! Common code for Picard assembly on the bulk mesh or on submeshes.
    [[nodiscard]] std::exception_ptr assembleOnBulkMeshOrOnSubmeshCommon(
        Assembly::CommonAssemblyData const& assembly_data,
        AssembledMatrixCache& mat_cache, double const t, double const dt,
        std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
        int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        std::vector<std::size_t> const* const sorted_element_subset)
    {
        std::exception_ptr exception = nullptr;

        try
        {
            auto const& loc_asms = derived().local_assemblers_;

            pvma_.assemble(
                loc_asms,
                assembly_data.activeElementIDsSorted(sorted_element_subset)
                    .get(),
                dof_tables, t, dt, x, x_prev, process_id, M, K, b);
        }
        catch (NumLib::AssemblyException const&)
        {
            exception = std::current_exception();
        }

        MathLib::LinAlg::finalizeAssembly(M);
        MathLib::LinAlg::finalizeAssembly(K);
        MathLib::LinAlg::finalizeAssembly(b);

        if (is_linear_)
        {
            DBUG("Saving global M, K, b for later reuse.");

            BaseLib::RunTime time_save;
            time_save.start();

            mat_cache.storeMKb(M, K, b);

            INFO("[time] Saving global M, K, b took {:g} s",
                 time_save.elapsed());
        }

        return exception;
    }

    std::exception_ptr assembleOnBulkMesh(
        double const t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
        int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        std::vector<std::size_t> const* const sorted_element_subset)
    {
        if (!bulk_mesh_matrix_cache_.hasMKb())
        {
            return assembleOnBulkMeshOrOnSubmeshCommon(
                *bulk_mesh_assembly_data_, bulk_mesh_matrix_cache_, t, dt, x,
                x_prev, dof_tables, process_id, M, K, b, sorted_element_subset);
        }

        DBUG("Reusing saved global M, K, b.");

        BaseLib::RunTime time_restore;
        time_restore.start();

        // restore bulk mesh matrices
        auto const [M_, K_, b_] = bulk_mesh_matrix_cache_.MKb();
        MathLib::LinAlg::copy(M_, M);
        MathLib::LinAlg::copy(K_, K);
        MathLib::LinAlg::copy(b_, b);

        INFO("[time] Restoring global M, K, b took {:g} s",
             time_restore.elapsed());

        return nullptr;
    }

    [[nodiscard]] std::exception_ptr assembleOnSubmeshes(
        double const t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
        int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        std::vector<std::size_t> const* const sorted_element_subset,
        bool const copy_residua_to_mesh)
    {
        auto const mat_spec = derived().getMatrixSpecifications(process_id);
        MKbFromCacheOrFromThis cached_or_this;

        std::exception_ptr exception = nullptr;
        for (auto const& [sad, mat_cache] :
             ranges::zip_view(submesh_assembly_data_ | ranges::views::all,
                              submesh_matrix_cache_ | ranges::views::all))
        {
            auto [from_where, M_submesh, K_submesh, b_submesh] =
                cached_or_this.MKb(mat_cache, mat_spec);

            if (from_where != MKbFromCacheOrFromThis::From::Cache)
            {
                // nothing cached, assembly necessary
                M_submesh.setZero();
                K_submesh.setZero();
                b_submesh.setZero();

                exception = assembleOnBulkMeshOrOnSubmeshCommon(
                    sad, mat_cache, t, dt, x, x_prev, dof_tables, process_id,
                    M_submesh, K_submesh, b_submesh, sorted_element_subset);
            }

            // TODO: This operates on the entire global vector b once
            // for each submesh. If the number of submeshes is high,
            // that might not perform well.
            MathLib::LinAlg::axpy(M, 1.0, M_submesh);
            MathLib::LinAlg::axpy(K, 1.0, K_submesh);
            MathLib::LinAlg::axpy(b, 1.0, b_submesh);

            if (copy_residua_to_mesh)
            {
                // Note: computeResiduum() uses bwd/fwd Euler scheme.
                GlobalVector const neg_res_submesh =
                    computeResiduum(dt, *x[process_id], *x_prev[process_id],
                                    M_submesh, K_submesh, b_submesh);

                AssemblyMixinBase::copyResiduumVectorsToSubmesh(
                    process_id, neg_res_submesh, *(dof_tables[process_id]),
                    sad);
            }
        }

        return exception;
    }

    //! Common code for the Newton-Raphson assembly on the bulk mesh or on
    //! submeshes.
    [[nodiscard]] std::exception_ptr
    assembleWithJacobianOnBulkMeshOrOnSubmeshCommon(
        Assembly::CommonAssemblyData const& assembly_data, double const t,
        double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
        int const process_id, GlobalVector& b, GlobalMatrix& Jac,
        std::vector<std::size_t> const* const sorted_element_subset)
    {
        auto const& loc_asms = derived().local_assemblers_;
        std::exception_ptr exception = nullptr;

        try
        {
            pvma_.assembleWithJacobian(
                loc_asms,
                assembly_data.activeElementIDsSorted(sorted_element_subset)
                    .get(),
                dof_tables, t, dt, x, x_prev, process_id, b, Jac);
        }
        catch (NumLib::AssemblyException const&)
        {
            exception = std::current_exception();
        }

        MathLib::LinAlg::finalizeAssembly(b);

        return exception;
    }

    [[nodiscard]] std::exception_ptr assembleWithJacobianOnSubmeshes(
        double const t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
        int const process_id, GlobalVector& b, GlobalMatrix& Jac,
        std::vector<std::size_t> const* const sorted_element_subset,
        bool const /*copy_residua_to_mesh*/)
    {
        std::exception_ptr exception = nullptr;

        auto const& mat_spec = derived().getMatrixSpecifications(process_id);
        auto b_submesh =
            MathLib::MatrixVectorTraits<GlobalVector>::newInstance(mat_spec);

        for (auto const& sad : submesh_assembly_data_)
        {
            b_submesh->setZero();

            exception = assembleWithJacobianOnBulkMeshOrOnSubmeshCommon(
                sad, t, dt, x, x_prev, dof_tables, process_id, *b_submesh, Jac,
                sorted_element_subset);

            MathLib::LinAlg::axpy(b, 1.0, *b_submesh);

            // TODO enabling this check needs a refactoring of postTimestep(),
            // computeSecondaryVariable() etc. hooks.
            if (/*copy_residua_to_mesh*/ true)
            {
                AssemblyMixinBase::copyResiduumVectorsToSubmesh(
                    process_id, *b_submesh, *(dof_tables[process_id]), sad);
            }
        }

        return exception;
    }

    void preOutputPicard(double const t,
                         double const dt,
                         std::vector<GlobalVector*> const& x,
                         std::vector<GlobalVector*> const& x_prev,
                         int const process_id)
    {
        auto const matrix_specification =
            derived().getMatrixSpecifications(process_id);

        auto M = MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(
            matrix_specification);
        auto K = MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(
            matrix_specification);
        auto b = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
            matrix_specification);

        M->setZero();
        K->setZero();
        b->setZero();

        assemble(t, dt, x, x_prev, process_id, *M, *K, *b, nullptr, true);
    }

    void preOutputNewton(double const /*t*/,
                         double const /*dt*/,
                         std::vector<GlobalVector*> const& /*x*/,
                         std::vector<GlobalVector*> const& /*x_prev*/,
                         int const /*process_id*/)
    {
        // TODO enabling this needs a refactoring of postTimestep(),
        // computeSecondaryVariable() etc. hooks.
#if 0
        auto const matrix_specification =
            derived().getMatrixSpecifications(process_id);

        auto Jac = MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(
            matrix_specification);
        auto b = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
            matrix_specification);

        Jac->setZero();
        b->setZero();

        assembleWithJacobian(t, dt, x, x_prev, process_id, *b, *Jac, nullptr,
                             true);
#endif
    }
};

}  // namespace ProcessLib
