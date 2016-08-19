/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_CONVERGENCECRITERIONPERCOMPONENT_H
#define NUMLIB_CONVERGENCECRITERIONPERCOMPONENT_H

#include "ConvergenceCriterion.h"

namespace MeshLib
{
class Mesh;
}  // namespace MeshLib

namespace NumLib
{
class LocalToGlobalIndexMap;

class ConvergenceCriterionPerComponent : public ConvergenceCriterion
{
public:
    virtual void setDOFTable(NumLib::LocalToGlobalIndexMap const& dof_table,
                             MeshLib::Mesh const& mesh) = 0;
};

}  // namespace NumLib

#endif  // NUMLIB_CONVERGENCECRITERIONPERCOMPONENT_H
