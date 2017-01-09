/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_VECTORMATRIXASSEMBLER_H_
#define PROCESSLIB_VECTORMATRIXASSEMBLER_H_

#include <vector>
#include "NumLib/NumericsConfig.h"
#include "AbstractJacobianAssembler.h"

namespace NumLib
{
class LocalToGlobalIndexMap;
}  // NumLib

namespace ProcessLib
{
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

    //! Assembles\c M, \c K, and \c b.
    //! \remark Jacobian is not assembled here, see assembleWithJacobian().
    void assemble(std::size_t const mesh_item_id,
                  LocalAssemblerInterface& local_assembler,
                  NumLib::LocalToGlobalIndexMap const& dof_table,
                  double const t, GlobalVector const& x, GlobalMatrix& M,
                  GlobalMatrix& K, GlobalVector& b);

    //! Assembles \c M, \c K, \c b, and the Jacobian \c Jac of the residual.
    //! \note The Jacobian must be assembled.
    void assembleWithJacobian(std::size_t const mesh_item_id,
                              LocalAssemblerInterface& local_assembler,
                              NumLib::LocalToGlobalIndexMap const& dof_table,
                              const double t, GlobalVector const& x,
                              GlobalVector const& xdot, const double dxdot_dx,
                              const double dx_dx, GlobalMatrix& M,
                              GlobalMatrix& K, GlobalVector& b,
                              GlobalMatrix& Jac);

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

}   // namespace ProcessLib

#endif  // PROCESSLIB_VECTORMATRIXASSEMBLER_H_
