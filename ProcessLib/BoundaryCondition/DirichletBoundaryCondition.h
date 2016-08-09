/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_DIRICHLETBC_H
#define PROCESS_LIB_DIRICHLETBC_H

#include "NumLib/IndexValueVector.h"
#include "BoundaryCondition.h"

namespace ProcessLib
{
class DirichletBoundaryCondition : public BoundaryCondition
{
public:
    virtual NumLib::IndexValueVector<GlobalIndexType> getBCValues() = 0;

    void apply(const double /*t*/,
               GlobalVector const& /*x*/,
               GlobalMatrix& /*K*/,
               GlobalVector& /*b*/) override final
    {
        // Do nothing. Dirichlet BCs are handled specially.
    }
};

}  // namespace ProcessLib

#endif  // PROCESS_LIB_DIRICHLETBC_H
