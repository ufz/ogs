/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ConvergenceCriterion.h"

namespace MeshLib
{
class Mesh;
}  // namespace MeshLib

namespace NumLib
{
class LocalToGlobalIndexMap;

//! Interface for applying a convergence criterion individually to each
//! component of a multi-component solution or residual vector.
//! Component here means sub-vector, not single scalar vector entry.
class ConvergenceCriterionPerComponent : public ConvergenceCriterion
{
public:
    //! Sets the d.o.f. table used to extract data for a specific component.
    virtual void setDOFTable(NumLib::LocalToGlobalIndexMap const& dof_table,
                             MeshLib::Mesh const& mesh) = 0;
};

}  // namespace NumLib
