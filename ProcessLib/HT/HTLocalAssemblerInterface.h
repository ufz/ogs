/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   HTLocalAssemblerInterface.h
 *  Created on October 11, 2017, 1:35 PM
 */

#pragma once

#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

namespace ProcessLib
{
struct CoupledSolutionsForStaggeredScheme;

namespace HT
{
template <typename NodalRowVectorType, typename GlobalDimNodalMatrixType>
struct IntegrationPointData final
{
    IntegrationPointData(NodalRowVectorType N_,
                         GlobalDimNodalMatrixType dNdx_,
                         double const& integration_weight_)
        : N(std::move(N_)),
          dNdx(std::move(dNdx_)),
          integration_weight(integration_weight_)
    {
    }

    NodalRowVectorType const N;
    GlobalDimNodalMatrixType const dNdx;
    double const integration_weight;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

class HTLocalAssemblerInterface : public ProcessLib::LocalAssemblerInterface,
                                  public NumLib::ExtrapolatableElement
{
public:
    HTLocalAssemblerInterface() : _coupled_solutions(nullptr) {}
    void setStaggeredCoupledSolutions(
        std::size_t const /*mesh_item_id*/,
        CoupledSolutionsForStaggeredScheme* const coupling_term)
    {
        _coupled_solutions = coupling_term;
    }

    virtual std::vector<double> const& getIntPtDarcyVelocity(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const = 0;

protected:
    // TODO: remove _coupled_solutions or move integration point data from local
    // assembler class to a new class to make local assembler unique for each
    // process.
    /** Pointer to CoupledSolutionsForStaggeredScheme that is set in a member of
     *  Process class, setCoupledTermForTheStaggeredSchemeToLocalAssemblers.
     *  It is used for calculate the secondary variables like velocity for
     *  coupled processes.
     */
    CoupledSolutionsForStaggeredScheme* _coupled_solutions;
};

}  // namespace HT
}  // namespace ProcessLib
