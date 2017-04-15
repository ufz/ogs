/**
 * \file
 * \author Karsten Rink
 * \date   2012-01-04
 * \brief  Implementation of the DirectConditionGenerator class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <fstream>
#include <memory>

#include <logog/include/logog.hpp>

#include "DirectConditionGenerator.h"

#include "Applications/FileIO/AsciiRasterInterface.h"

#include "Raster.h"
#include "MeshSurfaceExtraction.h"
#include "Mesh.h"
#include "MeshLib/Node.h"

#include <cmath>
#include <limits>

const std::vector<std::pair<std::size_t, double>>&
DirectConditionGenerator::directToSurfaceNodes(const MeshLib::Mesh& mesh,
                                               const std::string& filename)
{
    if (_direct_values.empty())
    {
        GeoLib::Raster* raster(
            FileIO::AsciiRasterInterface::readRaster(filename));
        if (!raster)
        {
            ERR("Error in DirectConditionGenerator::directToSurfaceNodes() - "
                "could not load raster file.");
            return _direct_values;
        }

        const MathLib::Vector3 dir(0, 0, -1);
        const std::vector<MeshLib::Node*> surface_nodes(
            MeshLib::MeshSurfaceExtraction::getSurfaceNodes(mesh, dir, 90));
        const double no_data(raster->getHeader().no_data);
        _direct_values.reserve(surface_nodes.size());
        for (auto const* surface_node : surface_nodes)
        {
            double val(raster->getValueAtPoint(*surface_node));
            val = (val == no_data) ? 0 : val;
            _direct_values.push_back(
                std::make_pair(surface_node->getID(), val));
        }
        delete raster;

        std::for_each(surface_nodes.begin(), surface_nodes.end(),
                      std::default_delete<MeshLib::Node>());
    }
    else
        ERR("Error in DirectConditionGenerator::directToSurfaceNodes() - Data "
            "vector contains outdated values.");

    return _direct_values;
}

const std::vector< std::pair<std::size_t,double> >& DirectConditionGenerator::directWithSurfaceIntegration(MeshLib::Mesh &mesh, const std::string &filename, double scaling)
{
    if (!_direct_values.empty()) {
        ERR("Error in DirectConditionGenerator::directWithSurfaceIntegration()"
            "- Data vector contains outdated values...");
        return _direct_values;
    }

    std::unique_ptr<GeoLib::Raster> raster(
        FileIO::AsciiRasterInterface::readRaster(filename));
    if (!raster) {
        ERR("Error in DirectConditionGenerator::directWithSurfaceIntegration()"
            "- could not load raster file.");
        return _direct_values;
    }

    MathLib::Vector3 const dir(0.0, 0.0, -1.0);
    double const angle(90);
    std::string const prop_name("OriginalSubsurfaceNodeIDs");
    std::unique_ptr<MeshLib::Mesh> surface_mesh(
        MeshLib::MeshSurfaceExtraction::getMeshSurface(
            mesh, dir, angle, prop_name));

    std::vector<double> node_area_vec =
        MeshLib::MeshSurfaceExtraction::getSurfaceAreaForNodes(*surface_mesh);
    const std::vector<MeshLib::Node*> &surface_nodes(surface_mesh->getNodes());
    const std::size_t nNodes(surface_mesh->getNumberOfNodes());
    const double no_data(raster->getHeader().no_data);

    auto const* const node_id_pv =
        surface_mesh->getProperties().getPropertyVector<int>(prop_name);
    if (!node_id_pv)
    {
        ERR(
            "Need subsurface node ids, but the property \"%s\" is not "
            "available.",
            prop_name.c_str());
        return _direct_values;
    }

    _direct_values.reserve(nNodes);
    for (std::size_t i=0; i<nNodes; ++i)
    {
        double val(raster->getValueAtPoint(*surface_nodes[i]));
        val = (val == no_data) ? 0 : ((val*node_area_vec[i])/scaling);
        _direct_values.emplace_back((*node_id_pv)[i], val);
    }

    return _direct_values;
}


int DirectConditionGenerator::writeToFile(const std::string &name) const
{
    std::ofstream out( name.c_str(), std::ios::out );

    if (out)
    {
        for (const auto& _direct_value : _direct_values)
            out << _direct_value.first << "\t" << _direct_value.second << "\n";

        out.close();
    }
    return 0;
}
