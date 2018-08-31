/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <logog/include/logog.hpp>

#include <memory>
#include <string>

#include "MeshLib/IO/readMeshFromFile.h"
// TODO used for output, if output classes are ready this has to be changed
#include "MeshLib/IO/writeMeshToFile.h"

#include "ProcessLib/CalculateSurfaceFlux/CalculateSurfaceFlux.h"
#include "ProcessLib/CalculateSurfaceFlux/ParseCalculateSurfaceFluxData.h"

namespace ProcessLib
{
struct Balance
{
    Balance(std::string&& balance_mesh_name,
            std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
            std::string&& balance_property_vector_name,
            std::string&& balance_output_mesh_file_name)
        : surface_mesh(*BaseLib::findElementOrError(
              meshes.begin(), meshes.end(),
              [&balance_mesh_name](auto const& m) {
                  return balance_mesh_name == m->getName();
              },
              "Expected to find a mesh named " + balance_mesh_name +
                  " for the balance calculation.")),
          mesh_name(std::move(balance_mesh_name)),
          property_vector_name(std::move(balance_property_vector_name)),
          output_mesh_file_name(std::move(balance_output_mesh_file_name))
    {
        DBUG(
            "read balance meta data:\n\tbalance mesh:\"%s\"\n\tproperty name: "
            "\"%s\"\n\toutput to: \"%s\"",
            mesh_name.c_str(), property_vector_name.c_str(),
            output_mesh_file_name.c_str());
    }

    static std::unique_ptr<Balance> createBalance(
        BaseLib::ConfigTree const& config,
        std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
        std::string const& output_directory)
    {
        std::string mesh_name;  // surface mesh the balance will computed on
        std::string balance_pv_name;
        std::string balance_out_fname;
        ProcessLib::parseCalculateSurfaceFluxData(
            config, mesh_name, balance_pv_name, balance_out_fname);

        if (mesh_name.empty())
        {
            return std::unique_ptr<ProcessLib::Balance>(nullptr);
        }
        balance_out_fname =
            BaseLib::copyPathToFileName(balance_out_fname, output_directory);
        return std::make_unique<Balance>(std::move(mesh_name), meshes,
                                         std::move(balance_pv_name),
                                         std::move(balance_out_fname));
    }

    void integrate(GlobalVector const& x, double const t, Process const& p,
                   int const process_id, int const integration_order,
                   MeshLib::Mesh const& bulk_mesh)
    {
        auto* const balance_pv = MeshLib::getOrCreateMeshProperty<double>(
            surface_mesh, property_vector_name, MeshLib::MeshItemType::Cell, 1);
        // initialise the PropertyVector pv with zero values
        std::fill(balance_pv->begin(), balance_pv->end(), 0.0);
        auto balance = ProcessLib::CalculateSurfaceFlux(
            surface_mesh,
            p.getProcessVariables(process_id)[0].get().getNumberOfComponents(),
            integration_order);

        balance.integrate(
            x, *balance_pv, t, bulk_mesh,
            [&p](std::size_t const element_id, MathLib::Point3d const& pnt,
                 double const t, GlobalVector const& x) {
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
