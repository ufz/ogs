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

#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "MeshLib/IO/readMeshFromFile.h"
// TODO used for output, if output classes are ready this has to be changed
#include "MeshLib/IO/writeMeshToFile.h"

#include "ProcessLib/SurfaceFlux/SurfaceFlux.h"

namespace ProcessLib
{
struct SurfaceFluxData
{
    SurfaceFluxData(
        std::string&& surfaceflux_mesh_name,
        std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
        std::string&& surfaceflux_property_vector_name,
        std::string&& surfaceflux_output_mesh_file_name)
        : surface_mesh(*BaseLib::findElementOrError(
              meshes.begin(), meshes.end(),
              [&surfaceflux_mesh_name](auto const& m) {
                  return surfaceflux_mesh_name == m->getName();
              },
              "Expected to find a mesh named " + surfaceflux_mesh_name +
                  " for the surfaceflux calculation.")),
          mesh_name(std::move(surfaceflux_mesh_name)),
          property_vector_name(std::move(surfaceflux_property_vector_name)),
          output_mesh_file_name(std::move(surfaceflux_output_mesh_file_name))
    {
        DBUG(
            "read surfaceflux meta data:\n\tsurfaceflux "
            "mesh:'{:s}'\n\tproperty "
            "name: '{:s}'\n\toutput to: '{:s}'",
            mesh_name, property_vector_name, output_mesh_file_name);
    }

    static std::unique_ptr<ProcessLib::SurfaceFluxData> createSurfaceFluxData(
        BaseLib::ConfigTree const& calculatesurfaceflux_config,
        std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
        std::string const& output_directory)
    {
        auto mesh_name =
            //! \ogs_file_param{prj__processes__process__calculatesurfaceflux__mesh}
            calculatesurfaceflux_config.getConfigParameter<std::string>("mesh");
        auto surfaceflux_pv_name =
            //! \ogs_file_param{prj__processes__process__calculatesurfaceflux__property_name}
            calculatesurfaceflux_config.getConfigParameter<std::string>(
                "property_name");
        auto surfaceflux_out_fname =
            //! \ogs_file_param{prj__processes__process__calculatesurfaceflux__output_mesh}
            calculatesurfaceflux_config.getConfigParameter<std::string>(
                "output_mesh");

        if (mesh_name.empty())
        {
            return nullptr;
        }
        surfaceflux_out_fname = BaseLib::copyPathToFileName(
            surfaceflux_out_fname, output_directory);
        return std::make_unique<SurfaceFluxData>(
            std::move(mesh_name), meshes, std::move(surfaceflux_pv_name),
            std::move(surfaceflux_out_fname));
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

    void save(double const t) const
    {
        // TODO (TomFischer) output, if output classes are ready this has to be
        // changed
        std::string const fname =
            BaseLib::dropFileExtension(output_mesh_file_name) + "_t_" +
            std::to_string(t) + ".vtu";
        MeshLib::IO::writeMeshToFile(surface_mesh, fname);
    }

private:
    MeshLib::Mesh& surface_mesh;
    std::string const mesh_name;
    std::string const property_vector_name;
    std::string const output_mesh_file_name;
};
}  // namespace ProcessLib
