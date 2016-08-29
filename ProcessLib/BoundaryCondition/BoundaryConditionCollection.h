/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_BOUNDARYCONDITIONCOLLECTION_H
#define PROCESSLIB_BOUNDARYCONDITIONCOLLECTION_H

#include "DirichletBoundaryCondition.h"
#include "ProcessLib/ProcessVariable.h"

namespace ProcessLib
{

class BoundaryConditionCollection final
{
public:
    BoundaryConditionCollection(
        std::vector<std::unique_ptr<ParameterBase>> const& parameters)
        : _parameters(parameters)
    {
    }

    void apply(const double t, GlobalVector const& x, GlobalMatrix& K,
               GlobalVector& b);

    std::vector<NumLib::IndexValueVector<GlobalIndexType>> const*
    getKnownSolutions(double const /*t*/) const
    {
        // TODO time-dependent Dirichlet BCs.
        return &_dirichlet_bcs;
    }

    void addBCsForProcessVariables(
        std::vector<std::reference_wrapper<ProcessVariable>> const&
            process_variables,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        unsigned const integration_order);

private:
    std::vector<NumLib::IndexValueVector<GlobalIndexType>> _dirichlet_bcs;
    std::vector<std::unique_ptr<BoundaryCondition>> _boundary_conditions;
    std::vector<std::unique_ptr<ParameterBase>> const& _parameters;
};


}  // ProcessLib

#endif  // PROCESSLIB_BOUNDARYCONDITIONCOLLECTION_H
