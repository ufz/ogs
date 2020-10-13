/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <string>

#include "BaseLib/Logging.h"
#include "ProcessLib/SurfaceFlux/SurfaceFlux.h"

namespace ProcessLib
{
struct SurfaceFluxData final
{
    SurfaceFluxData(MeshLib::Mesh& surfaceflux_mesh,
                    std::string&& surfaceflux_property_vector_name)
        : surface_mesh(surfaceflux_mesh),
          property_vector_name(std::move(surfaceflux_property_vector_name))
    {
    }

    static std::unique_ptr<ProcessLib::SurfaceFluxData> createSurfaceFluxData(
        BaseLib::ConfigTree const& calculatesurfaceflux_config,
        std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes)
    {
        auto mesh_name =
            //! \ogs_file_param{prj__processes__process__calculatesurfaceflux__mesh}
            calculatesurfaceflux_config.getConfigParameter<std::string>("mesh");
        auto surfaceflux_pv_name =
            //! \ogs_file_param{prj__processes__process__calculatesurfaceflux__property_name}
            calculatesurfaceflux_config.getConfigParameter<std::string>(
                "property_name");
        if (mesh_name.empty())
        {
            return nullptr;
        }

        DBUG(
            "Read surfaceflux meta data:\n\tmesh:'{:s}'\n\tproperty name: "
            "'{:s}'\n",
            mesh_name, surfaceflux_pv_name);

        auto& surface_mesh = *BaseLib::findElementOrError(
            meshes.begin(), meshes.end(),
            [&mesh_name](auto const& m) { return mesh_name == m->getName(); },
            "Expected to find a mesh named " + mesh_name +
                " for the surfaceflux calculation.");

        return std::make_unique<SurfaceFluxData>(
            surface_mesh, std::move(surfaceflux_pv_name));
    }

    void integrate(std::vector<GlobalVector*> const& x, double const t,
                   Process const& p, int const process_id,
                   int const integration_order, MeshLib::Mesh const& bulk_mesh,
                   std::vector<std::size_t> const& active_element_ids)
    {
        auto* const surfaceflux_pv = MeshLib::getOrCreateMeshProperty<double>(
            surface_mesh, property_vector_name, MeshLib::MeshItemType::Cell, 1);
        // initialise the PropertyVector pv with zero values
        std::fill(surfaceflux_pv->begin(), surfaceflux_pv->end(), 0.0);
        auto surfaceflux_process =
            ProcessLib::SurfaceFlux(surface_mesh,
                                    p.getProcessVariables(process_id)[0]
                                        .get()
                                        .getNumberOfGlobalComponents(),
                                    integration_order);

        surfaceflux_process.integrate(
            x, *surfaceflux_pv, t, bulk_mesh, active_element_ids,
            [&p](std::size_t const element_id, MathLib::Point3d const& pnt,
                 double const t, std::vector<GlobalVector*> const& x) {
                return p.getFlux(element_id, pnt, t, x);
            });
    }

private:
    MeshLib::Mesh& surface_mesh;
    std::string const property_vector_name;
};
}  // namespace ProcessLib
