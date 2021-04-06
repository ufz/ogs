/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshGenerator.h"

#include <memory>
#include <numeric>

#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Pyramid.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/MeshEditing/AddLayerToMesh.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
std::vector<MeshLib::Node*> MeshGenerator::generateRegularNodes(
    const std::vector<const std::vector<double>*>& vec_xyz_coords,
    const MathLib::Point3d& origin)
{
    auto const shift_coordinates = [](auto const& in, auto& out,
                                      auto const& shift) {
        std::transform(in.begin(), in.end(), std::back_inserter(out),
                       [&shift](auto const& v) { return v + shift; });
    };
    std::array<std::vector<double>, 3> coords;
    for (std::size_t i = 0; i < 3; ++i)
    {
        coords[i].reserve(vec_xyz_coords[i]->size());
        shift_coordinates(*vec_xyz_coords[i], coords[i], origin[i]);
    }

    std::vector<Node*> nodes;
    nodes.reserve(coords[0].size() * coords[1].size() * coords[2].size());

    for (auto const z : coords[2])
    {
        for (auto const y : coords[1])
        {
            std::transform(
                coords[0].begin(), coords[0].end(), std::back_inserter(nodes),
                [&y, &z](double const& x) { return new Node(x, y, z); });
        }
    }
    return nodes;
}

std::vector<MeshLib::Node*> MeshGenerator::generateRegularNodes(
    const std::vector<double>& vec_x_coords, const MathLib::Point3d& origin)
{
    std::vector<const std::vector<double>*> vec_xyz_coords;
    vec_xyz_coords.push_back(&vec_x_coords);
    std::vector<double> dummy(1, 0.0);
    for (unsigned i = vec_xyz_coords.size() - 1; i < 3u; i++)
    {
        vec_xyz_coords.push_back(&dummy);
    }
    return generateRegularNodes(vec_xyz_coords, origin);
}

std::vector<MeshLib::Node*> MeshGenerator::generateRegularNodes(
    std::vector<double>& vec_x_coords,
    std::vector<double>& vec_y_coords,
    const MathLib::Point3d& origin)
{
    std::vector<const std::vector<double>*> vec_xyz_coords;
    vec_xyz_coords.push_back(&vec_x_coords);
    vec_xyz_coords.push_back(&vec_y_coords);
    std::vector<double> dummy(1, 0.0);
    for (unsigned i = vec_xyz_coords.size() - 1; i < 3u; i++)
    {
        vec_xyz_coords.push_back(&dummy);
    }
    return generateRegularNodes(vec_xyz_coords, origin);
}

std::vector<MeshLib::Node*> MeshGenerator::generateRegularNodes(
    std::vector<double>& vec_x_coords,
    std::vector<double>& vec_y_coords,
    std::vector<double>& vec_z_coords,
    const MathLib::Point3d& origin)
{
    std::vector<const std::vector<double>*> vec_xyz_coords;
    vec_xyz_coords.push_back(&vec_x_coords);
    vec_xyz_coords.push_back(&vec_y_coords);
    vec_xyz_coords.push_back(&vec_z_coords);
    return generateRegularNodes(vec_xyz_coords, origin);
}

std::vector<MeshLib::Node*> MeshGenerator::generateRegularNodes(
    const std::array<unsigned, 3>& n_cells,
    const std::array<double, 3>& cell_size,
    const MathLib::Point3d& origin)
{
    std::vector<Node*> nodes;
    nodes.reserve((n_cells[0] + 1) * (n_cells[1] + 1) * (n_cells[2] + 1));

    for (std::size_t i = 0; i < n_cells[2] + 1; i++)
    {
        const double z(origin[2] + cell_size[2] * i);
        for (std::size_t j = 0; j < n_cells[1] + 1; j++)
        {
            const double y(origin[1] + cell_size[1] * j);
            for (std::size_t k = 0; k < n_cells[0] + 1; k++)
            {
                nodes.push_back(new Node(origin[0] + cell_size[0] * k, y, z));
            }
        }
    }
    return nodes;
}

namespace MeshGenerator
{
std::vector<MeshLib::Node*> generateRegularPyramidTopNodes(
    std::vector<double> const& x_coords,
    std::vector<double> const& y_coords,
    std::vector<double> const& z_coords,
    const MathLib::Point3d& origin)
{
    std::vector<Node*> nodes;
    nodes.reserve((x_coords.size() - 1) * (y_coords.size() - 1) *
                  (z_coords.size() - 1));

    auto const n_x_coords = x_coords.size() - 1;
    auto const n_y_coords = y_coords.size() - 1;
    auto const n_z_coords = z_coords.size() - 1;

    for (std::size_t i = 0; i < n_z_coords; i++)
    {
        const double z((z_coords[i] + z_coords[i + 1]) / 2 + origin[2]);
        for (std::size_t j = 0; j < n_y_coords; j++)
        {
            const double y((y_coords[j] + y_coords[j + 1]) / 2 + origin[1]);
            for (std::size_t k = 0; k < n_x_coords; k++)
            {
                const double x((x_coords[k] + x_coords[k + 1]) / 2 + origin[0]);
                nodes.push_back(new Node(x, y, z));
            }
        }
    }
    return nodes;
}
}  // end namespace MeshGenerator

Mesh* MeshGenerator::generateLineMesh(const double length,
                                      const std::size_t subdivision,
                                      const MathLib::Point3d& origin,
                                      std::string const& mesh_name)
{
    return generateLineMesh(subdivision, length / subdivision, origin,
                            mesh_name);
}

Mesh* MeshGenerator::generateLineMesh(const unsigned n_cells,
                                      const double cell_size,
                                      MathLib::Point3d const& origin,
                                      std::string const& mesh_name)
{
    return generateLineMesh(
        BaseLib::UniformSubdivision(n_cells * cell_size, n_cells), origin,
        mesh_name);
}

Mesh* MeshGenerator::generateLineMesh(const BaseLib::ISubdivision& div,
                                      MathLib::Point3d const& origin,
                                      std::string const& mesh_name)
{
    const std::vector<double> vec_x(div());
    std::vector<Node*> nodes(generateRegularNodes(vec_x, origin));

    // elements
    const std::size_t n_cells = nodes.size() - 1;
    std::vector<Element*> elements;
    elements.reserve(n_cells);

    for (std::size_t i = 0; i < n_cells; i++)
    {
        elements.push_back(new Line({nodes[i], nodes[i + 1]}));
    }

    return new Mesh(mesh_name, nodes, elements);
}

Mesh* MeshGenerator::generateRegularQuadMesh(const double length,
                                             const std::size_t subdivision,
                                             const MathLib::Point3d& origin,
                                             std::string const& mesh_name)
{
    return generateRegularQuadMesh(subdivision, subdivision,
                                   length / subdivision, length / subdivision,
                                   origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularQuadMesh(const double x_length,
                                             const double y_length,
                                             const std::size_t x_subdivision,
                                             const std::size_t y_subdivision,
                                             const MathLib::Point3d& origin,
                                             std::string const& mesh_name)
{
    return generateRegularQuadMesh(x_subdivision, y_subdivision,
                                   x_length / x_subdivision,
                                   y_length / y_subdivision, origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularQuadMesh(const unsigned n_x_cells,
                                             const unsigned n_y_cells,
                                             const double cell_size,
                                             MathLib::Point3d const& origin,
                                             std::string const& mesh_name)
{
    return generateRegularQuadMesh(n_x_cells, n_y_cells, cell_size, cell_size,
                                   origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularQuadMesh(const unsigned n_x_cells,
                                             const unsigned n_y_cells,
                                             const double cell_size_x,
                                             const double cell_size_y,
                                             MathLib::Point3d const& origin,
                                             std::string const& mesh_name)
{
    return generateRegularQuadMesh(
        BaseLib::UniformSubdivision(n_x_cells * cell_size_x, n_x_cells),
        BaseLib::UniformSubdivision(n_y_cells * cell_size_y, n_y_cells), origin,
        mesh_name);
}

Mesh* MeshGenerator::generateRegularQuadMesh(const BaseLib::ISubdivision& div_x,
                                             const BaseLib::ISubdivision& div_y,
                                             MathLib::Point3d const& origin,
                                             std::string const& mesh_name)
{
    std::vector<double> vec_x(div_x());
    std::vector<double> vec_y(div_y());
    std::vector<Node*> nodes(generateRegularNodes(vec_x, vec_y, origin));
    const unsigned n_x_nodes(vec_x.size());

    // elements
    std::vector<Element*> elements;
    const unsigned n_x_cells(vec_x.size() - 1);
    const unsigned n_y_cells(vec_y.size() - 1);
    elements.reserve(n_x_cells * n_y_cells);

    for (std::size_t j = 0; j < n_y_cells; j++)
    {
        const std::size_t offset_y1 = j * n_x_nodes;
        const std::size_t offset_y2 = (j + 1) * n_x_nodes;
        for (std::size_t k = 0; k < n_x_cells; k++)
        {
            elements.push_back(
                new Quad({nodes[offset_y1 + k], nodes[offset_y1 + k + 1],
                          nodes[offset_y2 + k + 1], nodes[offset_y2 + k]}));
        }
    }

    return new Mesh(mesh_name, nodes, elements);
}

Mesh* MeshGenerator::generateRegularHexMesh(const double length,
                                            const std::size_t subdivision,
                                            const MathLib::Point3d& origin,
                                            std::string const& mesh_name)
{
    return MeshGenerator::generateRegularHexMesh(
        subdivision, subdivision, subdivision, length / subdivision, origin,
        mesh_name);
}

Mesh* MeshGenerator::generateRegularHexMesh(const double x_length,
                                            const double y_length,
                                            const double z_length,
                                            const std::size_t x_subdivision,
                                            const std::size_t y_subdivision,
                                            const std::size_t z_subdivision,
                                            const MathLib::Point3d& origin,
                                            std::string const& mesh_name)
{
    return MeshGenerator::generateRegularHexMesh(
        x_subdivision, y_subdivision, z_subdivision, x_length / x_subdivision,
        y_length / y_subdivision, z_length / z_subdivision, origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularHexMesh(const unsigned n_x_cells,
                                            const unsigned n_y_cells,
                                            const unsigned n_z_cells,
                                            const double cell_size,
                                            MathLib::Point3d const& origin,
                                            std::string const& mesh_name)
{
    return MeshGenerator::generateRegularHexMesh(
        n_x_cells, n_y_cells, n_z_cells, cell_size, cell_size, cell_size,
        origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularHexMesh(const unsigned n_x_cells,
                                            const unsigned n_y_cells,
                                            const unsigned n_z_cells,
                                            const double cell_size_x,
                                            const double cell_size_y,
                                            const double cell_size_z,
                                            MathLib::Point3d const& origin,
                                            std::string const& mesh_name)
{
    return generateRegularHexMesh(
        BaseLib::UniformSubdivision(n_x_cells * cell_size_x, n_x_cells),
        BaseLib::UniformSubdivision(n_y_cells * cell_size_y, n_y_cells),
        BaseLib::UniformSubdivision(n_z_cells * cell_size_z, n_z_cells), origin,
        mesh_name);
}

Mesh* MeshGenerator::generateRegularHexMesh(const BaseLib::ISubdivision& div_x,
                                            const BaseLib::ISubdivision& div_y,
                                            const BaseLib::ISubdivision& div_z,
                                            MathLib::Point3d const& origin,
                                            std::string const& mesh_name)
{
    std::vector<double> vec_x(div_x());
    std::vector<double> vec_y(div_y());
    std::vector<double> vec_z(div_z());
    std::vector<Node*> nodes(generateRegularNodes(vec_x, vec_y, vec_z, origin));

    const unsigned n_x_nodes(vec_x.size());
    const unsigned n_y_nodes(vec_y.size());
    const unsigned n_x_cells(vec_x.size() - 1);
    const unsigned n_y_cells(vec_y.size() - 1);
    const unsigned n_z_cells(vec_z.size() - 1);

    // elements
    std::vector<Element*> elements;
    elements.reserve(n_x_cells * n_y_cells * n_z_cells);

    for (std::size_t i = 0; i < n_z_cells; i++)
    {
        const std::size_t offset_z1 = i * n_x_nodes * n_y_nodes;  // bottom
        const std::size_t offset_z2 = (i + 1) * n_x_nodes * n_y_nodes;  // top
        for (std::size_t j = 0; j < n_y_cells; j++)
        {
            const std::size_t offset_y1 = j * n_x_nodes;
            const std::size_t offset_y2 = (j + 1) * n_x_nodes;
            for (std::size_t k = 0; k < n_x_cells; k++)
            {
                elements.push_back(
                    new Hex({// bottom
                             nodes[offset_z1 + offset_y1 + k],
                             nodes[offset_z1 + offset_y1 + k + 1],
                             nodes[offset_z1 + offset_y2 + k + 1],
                             nodes[offset_z1 + offset_y2 + k],
                             // top
                             nodes[offset_z2 + offset_y1 + k],
                             nodes[offset_z2 + offset_y1 + k + 1],
                             nodes[offset_z2 + offset_y2 + k + 1],
                             nodes[offset_z2 + offset_y2 + k]}));
            }
        }
    }

    return new Mesh(mesh_name, nodes, elements);
}

Mesh* MeshGenerator::generateRegularPyramidMesh(
    const BaseLib::ISubdivision& div_x,
    const BaseLib::ISubdivision& div_y,
    const BaseLib::ISubdivision& div_z,
    MathLib::Point3d const& origin,
    std::string const& mesh_name)
{
    std::vector<double> vec_x(div_x());
    std::vector<double> vec_y(div_y());
    std::vector<double> vec_z(div_z());
    std::vector<Node*> nodes(generateRegularNodes(vec_x, vec_y, vec_z, origin));
    std::vector<Node*> const top_nodes(
        generateRegularPyramidTopNodes(vec_x, vec_y, vec_z, origin));

    nodes.insert(nodes.end(), top_nodes.begin(), top_nodes.end());

    const unsigned n_x_nodes(vec_x.size());
    const unsigned n_y_nodes(vec_y.size());
    const unsigned n_z_nodes(vec_z.size());
    const unsigned n_x_cells(vec_x.size() - 1);
    const unsigned n_y_cells(vec_y.size() - 1);
    const unsigned n_z_cells(vec_z.size() - 1);

    // elements
    std::vector<Element*> elements;
    auto const top_node_offset(n_x_nodes * n_y_nodes * n_z_nodes);
    elements.reserve(n_x_cells * n_y_cells * n_z_cells);

    for (std::size_t i = 0; i < n_z_cells; i++)
    {
        const std::size_t offset_z1 = i * n_x_nodes * n_y_nodes;  // bottom
        const std::size_t offset_z2 = (i + 1) * n_x_nodes * n_y_nodes;  // top
        for (std::size_t j = 0; j < n_y_cells; j++)
        {
            const std::size_t offset_y1 = j * n_x_nodes;
            const std::size_t offset_y2 = (j + 1) * n_x_nodes;
            for (std::size_t k = 0; k < n_x_cells; k++)
            {
                // generate 6 pyramids within the virtual hexahedron cell
                int const pyramid_top_index(i * n_x_cells * n_y_cells +
                                            j * n_x_cells + k +
                                            top_node_offset);
                elements.push_back(
                    new Pyramid{{// bottom 'hexahedron' face
                                 nodes[offset_z1 + offset_y1 + k],
                                 nodes[offset_z1 + offset_y1 + k + 1],
                                 nodes[offset_z1 + offset_y2 + k + 1],
                                 nodes[offset_z1 + offset_y2 + k],
                                 // top
                                 nodes[pyramid_top_index]}});
                elements.push_back(
                    new Pyramid{{// top 'hexahedron' face
                                 nodes[offset_z2 + offset_y1 + k + 1],
                                 nodes[offset_z2 + offset_y1 + k],
                                 nodes[offset_z2 + offset_y2 + k],
                                 nodes[offset_z2 + offset_y2 + k + 1],
                                 // top of pyramid directed towards the bottom
                                 nodes[pyramid_top_index]}});
                elements.push_back(
                    new Pyramid{{// right 'hexahedron' face
                                 nodes[offset_z1 + offset_y1 + k + 1],
                                 nodes[offset_z2 + offset_y1 + k + 1],
                                 nodes[offset_z2 + offset_y2 + k + 1],
                                 nodes[offset_z1 + offset_y2 + k + 1],
                                 // top of pyramid directed towards the bottom
                                 nodes[pyramid_top_index]}});
                elements.push_back(
                    new Pyramid{{// left 'hexahedron' face
                                 nodes[offset_z2 + offset_y1 + k],
                                 nodes[offset_z1 + offset_y1 + k],
                                 nodes[offset_z1 + offset_y2 + k],
                                 nodes[offset_z2 + offset_y2 + k],
                                 // top of pyramid directed towards the bottom
                                 nodes[pyramid_top_index]}});
                elements.push_back(
                    new Pyramid{{// front 'hexahedron' face
                                 nodes[offset_z2 + offset_y1 + k],
                                 nodes[offset_z2 + offset_y1 + k + 1],
                                 nodes[offset_z1 + offset_y1 + k + 1],
                                 nodes[offset_z1 + offset_y1 + k],
                                 // top of pyramid directed towards the bottom
                                 nodes[pyramid_top_index]}});
                elements.push_back(
                    new Pyramid{{// back 'hexahedron' face
                                 nodes[offset_z1 + offset_y2 + k],
                                 nodes[offset_z1 + offset_y2 + k + 1],
                                 nodes[offset_z2 + offset_y2 + k + 1],
                                 nodes[offset_z2 + offset_y2 + k],
                                 // top of pyramid directed towards the bottom
                                 nodes[pyramid_top_index]}});
            }
        }
    }
    return new Mesh(mesh_name, nodes, elements);
}

Mesh* MeshGenerator::generateRegularTriMesh(const double length,
                                            const std::size_t subdivision,
                                            const MathLib::Point3d& origin,
                                            std::string const& mesh_name)
{
    return generateRegularTriMesh(subdivision, subdivision,
                                  length / subdivision, origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularTriMesh(const double x_length,
                                            const double y_length,
                                            const std::size_t x_subdivision,
                                            const std::size_t y_subdivision,
                                            const MathLib::Point3d& origin,
                                            std::string const& mesh_name)
{
    return generateRegularTriMesh(x_subdivision, y_subdivision,
                                  x_length / x_subdivision,
                                  y_length / y_subdivision, origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularTriMesh(const unsigned n_x_cells,
                                            const unsigned n_y_cells,
                                            const double cell_size,
                                            MathLib::Point3d const& origin,
                                            std::string const& mesh_name)
{
    return generateRegularTriMesh(n_x_cells, n_y_cells, cell_size, cell_size,
                                  origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularTriMesh(const unsigned n_x_cells,
                                            const unsigned n_y_cells,
                                            const double cell_size_x,
                                            const double cell_size_y,
                                            MathLib::Point3d const& origin,
                                            std::string const& mesh_name)
{
    return generateRegularTriMesh(
        BaseLib::UniformSubdivision(n_x_cells * cell_size_x, n_x_cells),
        BaseLib::UniformSubdivision(n_y_cells * cell_size_y, n_y_cells), origin,
        mesh_name);
}

Mesh* MeshGenerator::generateRegularTriMesh(const BaseLib::ISubdivision& div_x,
                                            const BaseLib::ISubdivision& div_y,
                                            MathLib::Point3d const& origin,
                                            std::string const& mesh_name)
{
    std::vector<double> vec_x(div_x());
    std::vector<double> vec_y(div_y());
    std::vector<Node*> nodes(generateRegularNodes(vec_x, vec_y, origin));
    const unsigned n_x_nodes(vec_x.size());
    const unsigned n_x_cells(vec_x.size() - 1);
    const unsigned n_y_cells(vec_y.size() - 1);

    // elements
    std::vector<Element*> elements;
    elements.reserve(n_x_cells * n_y_cells * 2);

    for (std::size_t j = 0; j < n_y_cells; j++)
    {
        const std::size_t offset_y1 = j * n_x_nodes;
        const std::size_t offset_y2 = (j + 1) * n_x_nodes;
        for (std::size_t k = 0; k < n_x_cells; k++)
        {
            elements.push_back(
                new Tri({nodes[offset_y1 + k], nodes[offset_y2 + k + 1],
                         nodes[offset_y2 + k]}));

            elements.push_back(
                new Tri({nodes[offset_y1 + k], nodes[offset_y1 + k + 1],
                         nodes[offset_y2 + k + 1]}));
        }
    }

    return new Mesh(mesh_name, nodes, elements);
}

Mesh* MeshGenerator::generateRegularPrismMesh(const double x_length,
                                              const double y_length,
                                              const double z_length,
                                              const std::size_t x_subdivision,
                                              const std::size_t y_subdivision,
                                              const std::size_t z_subdivision,
                                              const MathLib::Point3d& origin,
                                              std::string const& mesh_name)
{
    return generateRegularPrismMesh(
        x_subdivision, y_subdivision, z_subdivision, x_length / x_subdivision,
        y_length / y_subdivision, z_length / z_subdivision, origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularPrismMesh(const unsigned n_x_cells,
                                              const unsigned n_y_cells,
                                              const unsigned n_z_cells,
                                              const double cell_size,
                                              MathLib::Point3d const& origin,
                                              std::string const& mesh_name)
{
    return generateRegularPrismMesh(n_x_cells, n_y_cells, n_z_cells, cell_size,
                                    cell_size, cell_size, origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularPrismMesh(const unsigned n_x_cells,
                                              const unsigned n_y_cells,
                                              const unsigned n_z_cells,
                                              const double cell_size_x,
                                              const double cell_size_y,
                                              const double cell_size_z,
                                              MathLib::Point3d const& origin,
                                              std::string const& mesh_name)
{
    std::unique_ptr<MeshLib::Mesh> mesh(generateRegularTriMesh(
        n_x_cells, n_y_cells, cell_size_x, cell_size_y, origin, mesh_name));
    std::size_t const n_tris(mesh->getNumberOfElements());
    bool const copy_material_ids = false;
    for (std::size_t i = 0; i < n_z_cells; ++i)
    {
        mesh.reset(MeshLib::addLayerToMesh(*mesh, cell_size_z, mesh_name, true,
                                           copy_material_ids));
    }
    std::vector<std::size_t> elem_ids(n_tris);
    std::iota(elem_ids.begin(), elem_ids.end(), 0);
    return MeshLib::removeElements(*mesh, elem_ids, mesh_name);
}

Mesh* MeshGenerator::generateRegularTetMesh(const double x_length,
                                            const double y_length,
                                            const double z_length,
                                            const std::size_t x_subdivision,
                                            const std::size_t y_subdivision,
                                            const std::size_t z_subdivision,
                                            const MathLib::Point3d& origin,
                                            std::string const& mesh_name)
{
    return generateRegularTetMesh(
        x_subdivision, y_subdivision, z_subdivision, x_length / x_subdivision,
        y_length / y_subdivision, z_length / z_subdivision, origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularTetMesh(const unsigned n_x_cells,
                                            const unsigned n_y_cells,
                                            const unsigned n_z_cells,
                                            const double cell_size_x,
                                            const double cell_size_y,
                                            const double cell_size_z,
                                            MathLib::Point3d const& origin,
                                            std::string const& mesh_name)
{
    return generateRegularTetMesh(
        BaseLib::UniformSubdivision(n_x_cells * cell_size_x, n_x_cells),
        BaseLib::UniformSubdivision(n_y_cells * cell_size_y, n_y_cells),
        BaseLib::UniformSubdivision(n_z_cells * cell_size_z, n_z_cells), origin,
        mesh_name);
}

Mesh* MeshGenerator::generateRegularTetMesh(const BaseLib::ISubdivision& div_x,
                                            const BaseLib::ISubdivision& div_y,
                                            const BaseLib::ISubdivision& div_z,
                                            MathLib::Point3d const& origin,
                                            std::string const& mesh_name)
{
    std::vector<double> vec_x(div_x());
    std::vector<double> vec_y(div_y());
    std::vector<double> vec_z(div_z());
    std::vector<Node*> nodes(generateRegularNodes(vec_x, vec_y, vec_z, origin));

    const unsigned n_x_nodes(vec_x.size());
    const unsigned n_y_nodes(vec_y.size());
    const unsigned n_x_cells(vec_x.size() - 1);
    const unsigned n_y_cells(vec_y.size() - 1);
    const unsigned n_z_cells(vec_z.size() - 1);

    // elements
    std::vector<Element*> elements;
    elements.reserve(n_x_cells * n_y_cells * n_z_cells * 6);

    for (std::size_t i = 0; i < n_z_cells; i++)
    {
        const std::size_t offset_z1 = i * n_x_nodes * n_y_nodes;  // bottom
        const std::size_t offset_z2 = (i + 1) * n_x_nodes * n_y_nodes;  // top
        for (std::size_t j = 0; j < n_y_cells; j++)
        {
            const std::size_t offset_y1 = j * n_x_nodes;
            const std::size_t offset_y2 = (j + 1) * n_x_nodes;
            for (std::size_t k = 0; k < n_x_cells; k++)
            {
                // tet 1
                elements.push_back(
                    new Tet({// bottom
                             nodes[offset_z1 + offset_y1 + k],
                             nodes[offset_z1 + offset_y2 + k + 1],
                             nodes[offset_z1 + offset_y2 + k],
                             // top
                             nodes[offset_z2 + offset_y1 + k]}));
                // tet 2
                elements.push_back(
                    new Tet({// bottom
                             nodes[offset_z1 + offset_y2 + k + 1],
                             nodes[offset_z1 + offset_y2 + k],
                             // top
                             nodes[offset_z2 + offset_y1 + k],
                             nodes[offset_z2 + offset_y2 + k + 1]}));
                // tet 3
                elements.push_back(
                    new Tet({// bottom
                             nodes[offset_z1 + offset_y2 + k],
                             // top
                             nodes[offset_z2 + offset_y1 + k],
                             nodes[offset_z2 + offset_y2 + k + 1],
                             nodes[offset_z2 + offset_y2 + k]}));
                // tet 4
                elements.push_back(
                    new Tet({// bottom
                             nodes[offset_z1 + offset_y1 + k],
                             nodes[offset_z1 + offset_y1 + k + 1],
                             nodes[offset_z1 + offset_y2 + k + 1],
                             // top
                             nodes[offset_z2 + offset_y1 + k + 1]}));
                // tet 5
                elements.push_back(
                    new Tet({// bottom
                             nodes[offset_z1 + offset_y1 + k],
                             nodes[offset_z1 + offset_y2 + k + 1],
                             // top
                             nodes[offset_z2 + offset_y1 + k],
                             nodes[offset_z2 + offset_y1 + k + 1]}));
                // tet 6
                elements.push_back(
                    new Tet({// bottom
                             nodes[offset_z1 + offset_y2 + k + 1],
                             // top
                             nodes[offset_z2 + offset_y1 + k],
                             nodes[offset_z2 + offset_y1 + k + 1],
                             nodes[offset_z2 + offset_y2 + k + 1]}));
            }
        }
    }

    return new Mesh(mesh_name, nodes, elements);
}

MeshLib::Mesh* MeshGenerator::createSurfaceMesh(
    std::string const& mesh_name, MathLib::Point3d const& ll,
    MathLib::Point3d const& ur, std::array<std::size_t, 2> const& n_steps,
    const std::function<double(double, double)>& f)
{
    std::array<double, 2> step_size{{(ur[0] - ll[0]) / (n_steps[0] - 1),
                                     (ur[1] - ll[1]) / (n_steps[1] - 1)}};

    std::vector<MeshLib::Node*> nodes;
    for (std::size_t j(0); j < n_steps[1]; ++j)
    {
        for (std::size_t i(0); i < n_steps[0]; ++i)
        {
            std::size_t const id = i + j * n_steps[1];
            double const x = ll[0] + i * step_size[0];
            double const y = ll[1] + j * step_size[1];

            nodes.push_back(new MeshLib::Node({x, y, f(x, y)}, id));
        }
    }

    std::vector<MeshLib::Element*> sfc_eles;
    for (std::size_t j(0); j < n_steps[1] - 1; ++j)
    {
        for (std::size_t i(0); i < n_steps[0] - 1; ++i)
        {
            std::size_t id_ll(i + j * n_steps[0]);
            std::size_t id_lr(i + 1 + j * n_steps[0]);
            std::size_t id_ul(i + (j + 1) * n_steps[0]);
            std::size_t id_ur(i + 1 + (j + 1) * n_steps[0]);
            sfc_eles.push_back(
                new MeshLib::Tri({nodes[id_ll], nodes[id_lr], nodes[id_ur]}));
            sfc_eles.push_back(
                new MeshLib::Tri({nodes[id_ll], nodes[id_ur], nodes[id_ul]}));
        }
    }

    return new MeshLib::Mesh(mesh_name, nodes, sfc_eles);
}

}  // namespace MeshLib
