/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Parameter/Parameter.h"

namespace ProcessLib
{
class SourceTerm
{
public:
    explicit SourceTerm(
        const NumLib::LocalToGlobalIndexMap& source_term_dof_table)
        : _source_term_dof_table(source_term_dof_table)
    {
    }

    virtual void integrate(const double t, GlobalVector& b) const = 0;

    virtual ~SourceTerm() = default;

protected:
    NumLib::LocalToGlobalIndexMap const& _source_term_dof_table;
};

}  // namespace ProcessLib
