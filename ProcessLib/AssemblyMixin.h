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

#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"
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
        std::vector<std::reference_wrapper<MeshLib::PropertyVector<double>>>&&
            residuum_vectors);

    MeshLib::PropertyVector<std::size_t> const& bulk_element_ids;
    MeshLib::PropertyVector<std::size_t> const& bulk_node_ids;
    std::vector<std::size_t> active_element_ids;
    std::vector<std::reference_wrapper<MeshLib::PropertyVector<double>>>
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
        const int process_id,
        MeshLib::Mesh& bulk_mesh,
        std::vector<std::reference_wrapper<MeshLib::Mesh>> const& submeshes,
        std::vector<std::string> const& residuum_names,
        std::vector<std::reference_wrapper<ProcessVariable>> const& pvs);

    void updateActiveElements(ProcessLib::ProcessVariable const& pv);

    static void copyResiduumVectorsToBulkMesh(
        GlobalVector const& rhs,
        NumLib::LocalToGlobalIndexMap const& local_to_global_index_map,
        std::vector<std::reference_wrapper<MeshLib::PropertyVector<double>>>
            residuum_vectors);

    static void copyResiduumVectorsToSubmesh(
        GlobalVector const& rhs,
        NumLib::LocalToGlobalIndexMap const& local_to_global_index_map,
        SubmeshAssemblyData const& sad);

private:
    void updateActiveElementsImpl(ProcessLib::ProcessVariable const& pv);

protected:
    std::vector<SubmeshAssemblyData> submesh_assembly_data_;
    std::vector<std::reference_wrapper<MeshLib::PropertyVector<double>>>
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
        const int process_id,
        std::vector<std::reference_wrapper<MeshLib::Mesh>> const& submeshes,
        std::vector<std::string> const& residuum_names)
    {
        AssemblyMixinBase::initializeAssemblyOnSubmeshes(
            process_id, derived().getMesh(), submeshes, residuum_names,
            derived().getProcessVariables(process_id));
    }

    void updateActiveElements(const int process_id)
    {
        // convention: process variable 0 governs where assembly takes place
        // (active element IDs)
        ProcessLib::ProcessVariable const& pv =
            derived().getProcessVariables(process_id)[0];

        AssemblyMixinBase::updateActiveElements(pv);
    }

    // cppcheck-suppress functionStatic
    void assemble(const double /*t*/, double const /*dt*/,
                  std::vector<GlobalVector*> const& /*x*/,
                  std::vector<GlobalVector*> const& /*x_prev*/,
                  int const /*process_id*/, GlobalMatrix& /*M*/,
                  GlobalMatrix& /*K*/, GlobalVector& /*b*/)
    {
        /*
        DBUG("AssemblyMixin assemble(t={}, dt={}, process_id={}).", t, dt,
             process_id);

        assembleGeneric(&Assembly::ParallelVectorMatrixAssembler::assemble, t,
        dt, x, x_prev, process_id, M, K, b);
        */
        OGS_FATAL("Not yet implemented.");
    }

    void assembleWithJacobian(const double t, double const dt,
                              std::vector<GlobalVector*> const& x,
                              std::vector<GlobalVector*> const& x_prev,
                              int const process_id, GlobalMatrix& M,
                              GlobalMatrix& K, GlobalVector& b,
                              GlobalMatrix& Jac)
    {
        DBUG("AssemblyMixin assembleWithJacobian(t={}, dt={}, process_id={}).",
             t, dt, process_id);

        assembleGeneric(
            &Assembly::ParallelVectorMatrixAssembler::assembleWithJacobian, t,
            dt, x, x_prev, process_id, M, K, b, Jac);
    }

private:
    Process& derived() { return static_cast<Process&>(*this); }
    Process const& derived() const
    {
        return static_cast<Process const&>(*this);
    }

    /// Generic assembly routine covering both the case with and without
    /// Jacobian assembly.
    template <typename Method, typename... Jac>
    void assembleGeneric(Method global_assembler_method, const double t,
                         double const dt, std::vector<GlobalVector*> const& x,
                         std::vector<GlobalVector*> const& x_prev,
                         int const process_id, GlobalMatrix& M, GlobalMatrix& K,
                         GlobalVector& b, Jac&... jac_or_not_jac)
    {
        // TODO why not getDOFTables(x.size()); ?
        std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const
            dof_tables{*derived()._local_to_global_index_map};

        auto const& loc_asms = derived().local_assemblers_;

        if (!submesh_assembly_data_.empty())
        {
            auto& b_submesh = NumLib::GlobalVectorProvider::provider.getVector(
                b, b_submesh_id_);

            for (auto const& sad : submesh_assembly_data_)
            {
                b_submesh.setZero();

                (pvma_.*global_assembler_method)(
                    loc_asms, sad.active_element_ids, dof_tables, t, dt, x,
                    x_prev, process_id, M, K, b_submesh, jac_or_not_jac...);

                MathLib::LinAlg::axpy(b, 1.0, b_submesh);

                AssemblyMixinBase::copyResiduumVectorsToSubmesh(
                    b_submesh, dof_tables.front().get(), sad);
            }

            NumLib::GlobalVectorProvider::provider.releaseVector(b_submesh);
        }
        else
        {
            // convention: process variable 0 governs where assembly takes
            // place (active element IDs)
            ProcessLib::ProcessVariable const& pv =
                derived().getProcessVariables(process_id)[0];

            (pvma_.*global_assembler_method)(
                loc_asms, pv.getActiveElementIDs(), dof_tables, t, dt, x,
                x_prev, process_id, M, K, b, jac_or_not_jac...);
        }

        AssemblyMixinBase::copyResiduumVectorsToBulkMesh(
            b, dof_tables.front().get(), residuum_vectors_bulk_);
    }
};

}  // namespace ProcessLib
