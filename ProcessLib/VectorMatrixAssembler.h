/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#include "NumLib/NumericsConfig.h"
#include "AbstractJacobianAssembler.h"
#include "CoupledSolutionsForStaggeredScheme.h"

namespace NumLib
{
class LocalToGlobalIndexMap;
}  // NumLib

namespace ProcessLib
{
struct CoupledSolutionsForStaggeredScheme;

class LocalAssemblerInterface;

//! Utility class used to assemble global matrices and vectors.
//!
//! The methods of this class get the global matrices and vectors as input and
//! pass only local data on to the local assemblers.
class VectorMatrixAssembler final
{
public:
    explicit VectorMatrixAssembler(
        std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler);

    void preAssemble(const std::size_t mesh_item_id,
                     LocalAssemblerInterface& local_assembler,
                     const NumLib::LocalToGlobalIndexMap& dof_table,
                     const double t, double const dt, const GlobalVector& x);

    //! Assembles\c M, \c K, and \c b.
    //! \remark Jacobian is not assembled here, see assembleWithJacobian().
    void assemble(std::size_t const mesh_item_id,
                  LocalAssemblerInterface& local_assembler,
                  std::vector<std::reference_wrapper<
                      NumLib::LocalToGlobalIndexMap>> const& dof_tables,
                  double const t, double const dt,
                  std::vector<GlobalVector*> const& x,
                  std::vector<GlobalVector*> const& xdot, int const process_id,
                  GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b);

    //! Assembles \c M, \c K, \c b, and the Jacobian \c Jac of the residual.
    //! \note The Jacobian must be assembled.
    void assembleWithJacobian(
        std::size_t const mesh_item_id,
        LocalAssemblerInterface& local_assembler,
        std::vector<
            std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
            dof_tables,
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& xdot, int const process_id,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac);

private:
    // temporary data only stored here in order to avoid frequent memory
    // reallocations.
    std::vector<double> _local_M_data;
    std::vector<double> _local_K_data;
    std::vector<double> _local_b_data;
    std::vector<double> _local_Jac_data;

    //! Used to assemble the Jacobian.
    std::unique_ptr<AbstractJacobianAssembler> _jacobian_assembler;
};

}  // namespace ProcessLib
