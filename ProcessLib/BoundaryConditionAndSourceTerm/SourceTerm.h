// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
class SourceTermBase
{
public:
    virtual void integrate(const double t, GlobalVector const& x,
                           GlobalVector& b, GlobalMatrix* jac) const = 0;

    virtual ~SourceTermBase() = default;
};

class SourceTerm : public SourceTermBase
{
public:
    explicit SourceTerm(
        std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table)
        : _source_term_dof_table{std::move(source_term_dof_table)}
    {
    }

protected:
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> const _source_term_dof_table;
};

}  // namespace ProcessLib
