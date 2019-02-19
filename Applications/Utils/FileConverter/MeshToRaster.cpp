/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <memory>
#include <string>

#include <tclap/CmdLine.h>

#include <Applications/ApplicationsLib/LogogSetup.h>

#include "BaseLib/BuildInfo.h"
#include "GeoLib/AABB.h"
#include "GeoLib/AnalyticalGeometry.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSearch/MeshElementGrid.h"
#include "MeshLib/Node.h"

/// Returns the index of the element the given node is located in when projected
/// onto a mesh.
std::size_t getProjectedElementIndex(
    std::vector<const MeshLib::Element*> const& elems,
    MeshLib::Node const& node)
{
    auto is_right_of = [&node](MeshLib::Node const& a, MeshLib::Node const& b) {
        return GeoLib::getOrientationFast(node, a, b) ==
               GeoLib::Orientation::CW;
    };

    std::size_t const n_elems = elems.size();
    for (std::size_t i = 0; i < n_elems; ++i)
    {
        if (elems[i]->getGeomType() == MeshLib::MeshElemType::LINE)
        {
            continue;
        }
        auto const& a = *elems[i]->getNode(0);
        auto const& b = *elems[i]->getNode(1);
        if (is_right_of(a, b))
        {
            continue;
        }
        auto const& c = *elems[i]->getNode(2);
        if (is_right_of(b, c))
        {
            continue;
        }
        if (elems[i]->getGeomType() == MeshLib::MeshElemType::TRIANGLE)
        {
            if (is_right_of(c, a))
            {
                continue;
            }
        }
        if (elems[i]->getGeomType() == MeshLib::MeshElemType::QUAD)
        {
            auto const& d = *elems[i]->getNode(3);
            if (is_right_of(c, d))
            {
                continue;
            }
            if (is_right_of(d, a))
            {
                continue;
            }
        }
        return elems[i]->getID();
    }
    return std::numeric_limits<size_t>::max();
}

/// Returns the z-coordinate of a point projected onto the plane defined by a
/// mesh element.
double getElevation(MeshLib::Element const& elem, MeshLib::Node const& node)
{
    MathLib::Vector3 const n =
        MeshLib::FaceRule::getSurfaceNormal(&elem).getNormalizedVector();
    MeshLib::Node const& orig = *elem.getNode(0);
    MathLib::Point3d const v{
        {node[0] - orig[0], node[1] - orig[1], node[2] - orig[2]}};
    double const dist = n[0] * v[0] + n[1] * v[1] + n[2] * v[2];
    return node[2] - dist * n[2];
}

int main(int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd(
        "Mesh to raster converter.\n"
        "Rasterises a 2D mesh, pixel values are set to the elevation of a "
        "regular grid superimposed on the mesh. If no mesh element is located "
        "beneath  a pixel it is set to NODATA.\n\n"
        "OpenGeoSys-6 software, version " +
            BaseLib::BuildInfo::git_describe +
            ".\n"
            "Copyright (c) 2012-2019, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', BaseLib::BuildInfo::git_describe);
    TCLAP::ValueArg<std::string> input_arg("i", "input-file",
                                           "Mesh input file (*.vtu, *.msh)",
                                           true, "", "string");
    cmd.add(input_arg);
    TCLAP::ValueArg<std::string> output_arg(
        "o", "output-file", "Raster output file (*.asc)", true, "", "string");
    cmd.add(output_arg);
    TCLAP::ValueArg<double> cell_arg("c", "cellsize",
                                     "edge length of raster cells in result",
                                     false, 1, "real");
    cmd.add(cell_arg);

    cmd.parse(argc, argv);

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

    double const cellsize =
        (cell_arg.isSet()) ? cell_arg.getValue() : mesh->getMinEdgeLength();
    INFO("Cellsize set to %f", cellsize);

    std::vector<MeshLib::Node*> const& nodes_vec(mesh->getNodes());
    GeoLib::AABB const bounding_box(nodes_vec.begin(), nodes_vec.end());
    MathLib::Point3d const& min(bounding_box.getMinPoint());
    MathLib::Point3d const& max(bounding_box.getMaxPoint());
    std::size_t const n_cols =
        static_cast<std::size_t>(std::ceil((max[0] - min[0]) / cellsize));
    std::size_t const n_rows =
        static_cast<std::size_t>(std::ceil((max[1] - min[1]) / cellsize));

    // raster header
    std::ofstream out(output_arg.getValue());
    out << "ncols         " << n_cols + 1 << "\n";
    out << "nrows         " << n_rows + 1 << "\n";
    out << std::fixed << "xllcorner     " << (min[0] - cellsize / 2.0) << "\n";
    out << std::fixed << "yllcorner     " << (min[1] - cellsize / 2.0) << "\n";
    out << std::fixed << "cellsize      " << cellsize << "\n";
    out << "NODATA_value  "
        << "-9999\n";
    INFO("Writing raster with %d x %d pixels.", n_cols, n_rows);

    MeshLib::MeshElementGrid const grid(*mesh);
    double const max_edge(mesh->getMaxEdgeLength() + cellsize);

    // pixel values
    double x = min[0];
    double y = max[1];
    double const half_cell = cellsize / 2.0;
    for (std::size_t j = 0; j <= n_rows; ++j)
    {
        for (std::size_t i = 0; i <= n_cols; ++i)
        {
            MeshLib::Node const node(x, y, 0);
            MathLib::Point3d min_vol{{x - max_edge, y - max_edge,
                                      -std::numeric_limits<double>::max()}};
            MathLib::Point3d max_vol{{x + max_edge, y + max_edge,
                                      std::numeric_limits<double>::max()}};
            std::vector<const MeshLib::Element*> const elems =
                grid.getElementsInVolume(min_vol, max_vol);
            std::size_t const elem_idx = getProjectedElementIndex(elems, node);
            // centre of the pixel is located within a mesh element
            if (elem_idx != std::numeric_limits<std::size_t>::max())
                out << getElevation(*(mesh->getElement(elem_idx)), node) << " ";
            // centre of the pixel isn't located within a mesh element
            else
            {
                std::array<double, 4> const x_off{
                    {-half_cell, half_cell, -half_cell, half_cell}};
                std::array<double, 4> const y_off{
                    {-half_cell, -half_cell, half_cell, half_cell}};
                double sum(0);
                std::size_t nonzero_count(0);
                // test for all of the pixel's corner if any are within an
                // element
                for (std::size_t i = 0; i < 4; ++i)
                {
                    MeshLib::Node const node(x + x_off[i], y + y_off[i], 0);
                    std::size_t const corner_idx =
                        getProjectedElementIndex(elems, node);
                    if (corner_idx != std::numeric_limits<std::size_t>::max())
                    {
                        sum +=
                            getElevation(*(mesh->getElement(corner_idx)), node);
                        nonzero_count++;
                    }
                }
                // calculate pixel value as average of corner values
                if (nonzero_count > 0)
                    out << sum / nonzero_count << " ";
                // if none of the corners give a value, set pixel to NODATA
                else
                    out << "-9999 ";
            }
            x += cellsize;
        }
        x = min[0];
        y -= cellsize;
        out << "\n";
    }
    out.close();
    INFO("Result written to %s", output_arg.getValue().c_str());
    return 0;
}
