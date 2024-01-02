/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "BaseLib/ConfigTree-fwd.h"
#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace MeshLib
{
class Mesh;
}  // namespace MeshLib

namespace ProcessLib
{
class Process;

struct SurfaceFluxData final
{
    SurfaceFluxData(MeshLib::Mesh& surfaceflux_mesh,
                    std::string&& surfaceflux_property_vector_name);

    static std::unique_ptr<ProcessLib::SurfaceFluxData> createSurfaceFluxData(
        BaseLib::ConfigTree const& calculatesurfaceflux_config,
        std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes);

    void integrate(std::vector<GlobalVector*> const& x, double const t,
                   Process const& p, int const process_id,
                   int const integration_order, MeshLib::Mesh const& bulk_mesh,
                   std::vector<std::size_t> const& active_element_ids);

private:
    MeshLib::Mesh& surface_mesh;
    std::string const property_vector_name;
};
}  // namespace ProcessLib
