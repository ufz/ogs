/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
class SourceTerm
{
public:
    explicit SourceTerm(
        std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table)
        : _source_term_dof_table{std::move(source_term_dof_table)}
    {
    }

    virtual void integrate(const double t, GlobalVector const& x,
                           GlobalVector& b, GlobalMatrix* jac) const = 0;

    virtual ~SourceTerm() = default;

protected:
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> const _source_term_dof_table;
};

}  // namespace ProcessLib
