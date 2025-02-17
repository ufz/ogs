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

#include "BaseLib/MPI.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"
#include "NumLib/Exceptions.h"
#include "ProcessLib/Assembly/ParallelVectorMatrixAssembler.h"
#include "ProcessLib/ProcessVariable.h"
#include "VectorMatrixAssembler.h"

namespace ProcessLib
{
//! Data necessary for global equation system assembly on submeshes of the bulk
//! mesh.
struct SubmeshAssemblyData
{
    explicit SubmeshAssemblyData(
        MeshLib::Mesh const& mesh,
        std::vector<std::vector<
            std::reference_wrapper<MeshLib::PropertyVector<double>>>>&&
            residuum_vectors);

    MeshLib::PropertyVector<std::size_t> const& bulk_element_ids;
    MeshLib::PropertyVector<std::size_t> const& bulk_node_ids;
    std::vector<std::size_t> active_element_ids;
    std::vector<
        std::vector<std::reference_wrapper<MeshLib::PropertyVector<double>>>>
        residuum_vectors;
};

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
    explicit AssemblyMixinBase(AbstractJacobianAssembler& jacobian_assembler)
        : pvma_{jacobian_assembler} {};

    void initializeAssemblyOnSubmeshes(
        MeshLib::Mesh& bulk_mesh,
        std::vector<std::reference_wrapper<MeshLib::Mesh>> const& submeshes,
        std::vector<std::vector<std::string>> const& residuum_names,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>> const&
            pvs);

    void updateActiveElements(ProcessLib::Process const& process);

    static void copyResiduumVectorsToBulkMesh(
        GlobalVector const& rhs,
        NumLib::LocalToGlobalIndexMap const& local_to_global_index_map,
        std::vector<std::reference_wrapper<MeshLib::PropertyVector<double>>>
            residuum_vectors);

    static void copyResiduumVectorsToSubmesh(
        int const process_id,
        GlobalVector const& rhs,
        NumLib::LocalToGlobalIndexMap const& local_to_global_index_map,
        SubmeshAssemblyData const& sad);

private:
    void updateActiveElementsImpl(Process const& process);

protected:
    std::vector<SubmeshAssemblyData> submesh_assembly_data_;
    std::vector<
        std::vector<std::reference_wrapper<MeshLib::PropertyVector<double>>>>
        residuum_vectors_bulk_;

    /// ID of the b vector on submeshes, cf. NumLib::VectorProvider.
    std::size_t b_submesh_id_ = 0;

    Assembly::ParallelVectorMatrixAssembler pvma_;

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

    // cppcheck-suppress functionStatic
    void assemble(double const t, double const dt,
                  std::vector<GlobalVector*> const& /*x*/,
                  std::vector<GlobalVector*> const& /*x_prev*/,
                  int const process_id, GlobalMatrix& /*M*/,
                  GlobalMatrix& /*K*/, GlobalVector& /*b*/)
    {
        DBUG("AssemblyMixin assemble(t={}, dt={}, process_id={}).", t, dt,
             process_id);

        /// Implementation similar to assembleWithJacobian calling
        /// Assembly::ParallelVectorMatrixAssembler::assemble function.
        /// Residuum must be correctly computed.
        OGS_FATAL("AssemblyMixin for Picard scheme is not yet implemented.");
    }

    void assembleWithJacobian(double const t, double const dt,
                              std::vector<GlobalVector*> const& x,
                              std::vector<GlobalVector*> const& x_prev,
                              int const process_id, GlobalVector& b,
                              GlobalMatrix& Jac)
    {
        DBUG("AssemblyMixin assembleWithJacobian(t={}, dt={}, process_id={}).",
             t, dt, process_id);

        std::vector<NumLib::LocalToGlobalIndexMap const*> const dof_tables =
            derived().getDOFTables(x.size());

        auto const& loc_asms = derived().local_assemblers_;

        std::exception_ptr exception = nullptr;
        if (!submesh_assembly_data_.empty())
        {
            auto& b_submesh = NumLib::GlobalVectorProvider::provider.getVector(
                b, b_submesh_id_);

            for (auto const& sad : submesh_assembly_data_)
            {
                b_submesh.setZero();

                try
                {
                    pvma_.assembleWithJacobian(loc_asms, sad.active_element_ids,
                                               dof_tables, t, dt, x, x_prev,
                                               process_id, b_submesh, Jac);
                }
                catch (NumLib::AssemblyException const&)
                {
                    exception = std::current_exception();
                }

                MathLib::LinAlg::axpy(b, 1.0, b_submesh);

                AssemblyMixinBase::copyResiduumVectorsToSubmesh(
                    process_id, b_submesh, *(dof_tables[process_id]), sad);
            }

            NumLib::GlobalVectorProvider::provider.releaseVector(b_submesh);
        }
        else
        {
            try
            {
                pvma_.assembleWithJacobian(
                    loc_asms, derived().getActiveElementIDs(), dof_tables, t,
                    dt, x, x_prev, process_id, b, Jac);
            }
            catch (NumLib::AssemblyException const&)
            {
                exception = std::current_exception();
            }
        }

        MathLib::LinAlg::finalizeAssembly(b);
        MathLib::LinAlg::finalizeAssembly(Jac);

        if (BaseLib::MPI::anyOf(exception != nullptr))
        {
            if (exception)  // Only the rank with the exception rethrows...
            {
                std::rethrow_exception(exception);
            }
            // ... but all ranks quit.
            return;
        }

        AssemblyMixinBase::copyResiduumVectorsToBulkMesh(
            b, *(dof_tables[process_id]), residuum_vectors_bulk_[process_id]);
    }

private:
    Process& derived() { return static_cast<Process&>(*this); }
    Process const& derived() const
    {
        return static_cast<Process const&>(*this);
    }
};

}  // namespace ProcessLib
