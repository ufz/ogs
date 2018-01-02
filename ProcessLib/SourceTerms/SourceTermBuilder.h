/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NodalSourceTerm.h"

namespace GeoLib
{
class GeoObject;
}

namespace MeshLib
{
class Element;
class Mesh;
}

namespace MeshGeoToolsLib
{
class BoundaryElementsSearcher;
}

namespace NumLib
{
class LocalToGlobalIndexMap;
template <typename>
struct IndexValueVector;
}

namespace ProcessLib
{
struct SourceTermConfig;

class SourceTermBuilder
{
public:
    virtual ~SourceTermBuilder() = default;

    virtual std::unique_ptr<NodalSourceTerm> createSourceTerm(
        const SourceTermConfig& config,
        const NumLib::LocalToGlobalIndexMap& dof_table,
        const MeshLib::Mesh& mesh, const int variable_id,
        const unsigned integration_order, const unsigned shapefunction_order);

protected:
    virtual std::unique_ptr<NodalSourceTerm> createNodalSourceTerm(
        const SourceTermConfig& config,
        const NumLib::LocalToGlobalIndexMap& dof_table,
        const MeshLib::Mesh& mesh, const int variable_id,
        const unsigned integration_order, const unsigned shapefunction_order);
};

}  // ProcessLib
