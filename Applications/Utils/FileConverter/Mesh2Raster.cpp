/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include <filesystem>
#include <fstream>
#include <memory>
#include <string>

#include "BaseLib/MPI.h"
#include "GeoLib/AABB.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSearch/MeshElementGrid.h"
#include "MeshLib/Node.h"
#include "MeshToolsLib/MeshEditing/ProjectPointOnMesh.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Mesh to raster converter.\n"
        "Rasterises a 2D mesh, pixel values are set to the elevation of a "
        "regular grid superimposed on the mesh. If no mesh element is located "
        "beneath  a pixel it is set to NODATA.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<double> cell_arg("c", "cellsize",
                                     "edge length of raster cells in result",
                                     false, 1, "real");
    cmd.add(cell_arg);
    TCLAP::ValueArg<std::string> output_arg(
        "o", "output-file", "Raster output file (*.asc)", true, "", "string");
    cmd.add(output_arg);
    TCLAP::ValueArg<std::string> input_arg("i", "input-file",
                                           "Mesh input file (*.vtu, *.msh)",
                                           true, "", "string");
    cmd.add(input_arg);
    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);

    INFO("Rasterising mesh...");
    std::unique_ptr<MeshLib::Mesh> const mesh(
        MeshLib::IO::readMeshFromFile(input_arg.getValue()));
    if (mesh == nullptr)
    {
        ERR("Error reading mesh file.");
        return 1;
    }
    if (mesh->getDimension() != 2)
    {
        ERR("The programme requires a mesh containing two-dimensional elements "
            "(i.e. triangles or quadrilaterals.");
        return 2;
    }

    auto const edgeLengths = minMaxEdgeLength(mesh->getElements());
    double const cellsize =
        (cell_arg.isSet()) ? cell_arg.getValue() : edgeLengths.first;
    INFO("Cellsize set to {:f}", cellsize);

    double const max_edge = edgeLengths.second + cellsize;

    std::vector<MeshLib::Node*> const& nodes_vec(mesh->getNodes());
    GeoLib::AABB const bounding_box(nodes_vec.begin(), nodes_vec.end());
    auto const& [min, max] = bounding_box.getMinMaxPoints();
    auto const n_cols =
        static_cast<std::size_t>(std::ceil((max[0] - min[0]) / cellsize));
    auto const n_rows =
        static_cast<std::size_t>(std::ceil((max[1] - min[1]) / cellsize));
    double const half_cell = cellsize / 2.0;

    // raster header
    std::string output_name = output_arg.getValue();
    std::string const ext = ".asc";
    if (std::filesystem::path(output_name).extension() != ext)
    {
        WARN("Adding extension '*.asc' to output file name '{:s}'.",
             output_name);
        output_name += ext;
    }
    std::ofstream out(output_name);
    out << "ncols         " << n_cols << "\n";
    out << "nrows         " << n_rows << "\n";
    out << std::fixed << "xllcorner     " << (min[0] - half_cell) << "\n";
    out << std::fixed << "yllcorner     " << (min[1] - half_cell) << "\n";
    out << std::fixed << "cellsize      " << cellsize << "\n";
    out << "NODATA_value  "
        << "-9999\n";
    INFO("Writing raster with {:d} x {:d} pixels.", n_cols, n_rows);

    MeshLib::MeshElementGrid const grid(*mesh);

    for (std::size_t row = 0; row < n_rows; ++row)
    {
        double const y = max[1] - row * cellsize;
        for (std::size_t column = 0; column < n_cols; ++column)
        {
            // pixel values
            double const x = min[0] + column * cellsize;
            MeshLib::Node const node(x, y, 0);
            MathLib::Point3d min_vol{{x - max_edge, y - max_edge,
                                      -std::numeric_limits<double>::max()}};
            MathLib::Point3d max_vol{{x + max_edge, y + max_edge,
                                      std::numeric_limits<double>::max()}};
            std::vector<const MeshLib::Element*> const& elems =
                grid.getElementsInVolume(min_vol, max_vol);
            auto const* element =
                MeshToolsLib::ProjectPointOnMesh::getProjectedElement(elems,
                                                                      node);
            // centre of the pixel is located within a mesh element
            if (element != nullptr)
            {
                out << MeshToolsLib::ProjectPointOnMesh::getElevation(*element,
                                                                      node)
                    << " ";
            }
            else
            {
                std::array<double, 4> const x_off{
                    {-half_cell, half_cell, -half_cell, half_cell}};
                std::array<double, 4> const y_off{
                    {-half_cell, -half_cell, half_cell, half_cell}};
                double sum(0);
                std::size_t nonzero_count(0);
                // test all of the pixel's corners if there are any within an
                // element
                for (std::size_t i = 0; i < 4; ++i)
                {
                    MeshLib::Node const corner_node(x + x_off[i], y + y_off[i],
                                                    0);
                    auto const* corner_element =
                        MeshToolsLib::ProjectPointOnMesh::getProjectedElement(
                            elems, corner_node);
                    if (corner_element != nullptr)
                    {
                        sum += MeshToolsLib::ProjectPointOnMesh::getElevation(
                            *corner_element, node);
                        nonzero_count++;
                    }
                }
                if (nonzero_count > 0)
                {
                    // calculate pixel value as average of corner values
                    out << sum / nonzero_count << " ";
                }
                else
                {
                    // if none of the corners give a value, set pixel to NODATA
                    out << "-9999 ";
                }
            }
        }
        out << "\n";
    }
    out.close();
    INFO("Result written to {:s}", output_name);
    return 0;
}
