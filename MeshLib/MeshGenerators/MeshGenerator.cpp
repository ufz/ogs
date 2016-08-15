/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshGenerator.h"

#include <memory>
#include <numeric>

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/MeshEditing/AddLayerToMesh.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"

namespace MeshLib
{

std::vector<MeshLib::Node*> MeshGenerator::generateRegularNodes(
    const std::vector<const std::vector<double>*> &vec_xyz_coords,
    const MathLib::Point3d& origin)
{
    std::vector<Node*> nodes;
    nodes.reserve(vec_xyz_coords[0]->size()*vec_xyz_coords[1]->size()*vec_xyz_coords[2]->size());

    for (std::size_t i = 0; i < vec_xyz_coords[2]->size(); i++)
    {
        const double z ((*vec_xyz_coords[2])[i]+origin[2]);
        for (std::size_t j = 0; j < vec_xyz_coords[1]->size(); j++)
        {
            const double y ((*vec_xyz_coords[1])[j]+origin[1]);
            for (double const x : *vec_xyz_coords[0])
            {
                nodes.push_back (new Node(x+origin[0], y, z));
            }
        }
    }
    return nodes;
}

std::vector<MeshLib::Node*> MeshGenerator::generateRegularNodes(
    const std::vector<double> &vec_x_coords,
    const MathLib::Point3d& origin)
{
    std::vector<const std::vector<double>*> vec_xyz_coords;
    vec_xyz_coords.push_back(&vec_x_coords);
    std::vector<double> dummy(1,0.0);
    for (unsigned i=vec_xyz_coords.size()-1; i<3u; i++)
        vec_xyz_coords.push_back(&dummy);
    return generateRegularNodes(vec_xyz_coords, origin);
}

std::vector<MeshLib::Node*> MeshGenerator::generateRegularNodes(
    std::vector<double> &vec_x_coords,
    std::vector<double> &vec_y_coords,
    const MathLib::Point3d& origin)
{
    std::vector<const std::vector<double>*> vec_xyz_coords;
    vec_xyz_coords.push_back(&vec_x_coords);
    vec_xyz_coords.push_back(&vec_y_coords);
    std::vector<double> dummy(1,0.0);
    for (unsigned i=vec_xyz_coords.size()-1; i<3u; i++)
        vec_xyz_coords.push_back(&dummy);
    return generateRegularNodes(vec_xyz_coords, origin);
}

std::vector<MeshLib::Node*> MeshGenerator::generateRegularNodes(
    std::vector<double> &vec_x_coords,
    std::vector<double> &vec_y_coords,
    std::vector<double> &vec_z_coords,
    const MathLib::Point3d& origin)
{
    std::vector<const std::vector<double>*> vec_xyz_coords;
    vec_xyz_coords.push_back(&vec_x_coords);
    vec_xyz_coords.push_back(&vec_y_coords);
    vec_xyz_coords.push_back(&vec_z_coords);
    return generateRegularNodes(vec_xyz_coords, origin);
}

std::vector<MeshLib::Node*> MeshGenerator::generateRegularNodes(
    const std::array<unsigned,3> &n_cells,
    const std::array<double,3> &cell_size,
    const MathLib::Point3d& origin)
{
    std::vector<Node*> nodes;
    nodes.reserve((n_cells[0]+1)*(n_cells[1]+1)*(n_cells[2]+1));

    for (std::size_t i = 0; i < n_cells[2]+1; i++)
    {
        const double z (origin[2] + cell_size[2] * i);
        for (std::size_t j = 0; j < n_cells[1]+1; j++)
        {
            const double y (origin[1] + cell_size[1] * j);
            for (std::size_t k = 0; k < n_cells[0]+1; k++)
            {
                nodes.push_back (new Node(origin[0] + cell_size[0] * k, y, z));
            }
        }
    }
    return nodes;
}

Mesh* MeshGenerator::generateLineMesh(
    const double length,
    const std::size_t subdivision,
    const MathLib::Point3d& origin,
    std::string   const& mesh_name)
{
    return generateLineMesh(subdivision, length/subdivision, origin, mesh_name);
}

Mesh* MeshGenerator::generateLineMesh(
    const unsigned n_cells,
    const double   cell_size,
    MathLib::Point3d const& origin,
    std::string   const& mesh_name)
{
    return generateLineMesh(BaseLib::UniformSubdivision(n_cells*cell_size, n_cells), origin, mesh_name);
}

Mesh* MeshGenerator::generateLineMesh(
    const BaseLib::ISubdivision &div,
    MathLib::Point3d const& origin,
    std::string const& mesh_name)
{
    const std::vector<double> vec_x(div());
    std::vector<Node*> nodes(generateRegularNodes(vec_x, origin));

    //elements
    const std::size_t n_cells = nodes.size()-1;
    std::vector<Element*> elements;
    elements.reserve(n_cells);

    for (std::size_t i = 0; i < n_cells; i++)
    {
        std::array<Node*, 2> element_nodes;
        element_nodes[0] = nodes[i];
        element_nodes[1] = nodes[i + 1];
        elements.push_back (new Line(element_nodes));
    }

    return new Mesh(mesh_name, nodes, elements);
}

Mesh* MeshGenerator::generateRegularQuadMesh(
    const double length,
    const std::size_t subdivision,
    const MathLib::Point3d& origin,
    std::string const& mesh_name)
{
    return generateRegularQuadMesh(subdivision, subdivision,
        length/subdivision, length/subdivision, origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularQuadMesh(
    const double x_length,
    const double y_length,
    const std::size_t x_subdivision,
    const std::size_t y_subdivision,
    const MathLib::Point3d& origin,
    std::string const& mesh_name)
{
    return generateRegularQuadMesh(x_subdivision, y_subdivision,
        x_length/x_subdivision, y_length/y_subdivision, origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularQuadMesh(
    const unsigned n_x_cells,
    const unsigned n_y_cells,
    const double cell_size,
    MathLib::Point3d const& origin,
    std::string const& mesh_name)
{
    return generateRegularQuadMesh(n_x_cells, n_y_cells, cell_size, cell_size, origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularQuadMesh(
    const unsigned n_x_cells,
    const unsigned n_y_cells,
    const double cell_size_x,
    const double cell_size_y,
    MathLib::Point3d const& origin,
    std::string const& mesh_name)
{
    return generateRegularQuadMesh(BaseLib::UniformSubdivision(n_x_cells*cell_size_x, n_x_cells),
        BaseLib::UniformSubdivision(n_y_cells*cell_size_y, n_y_cells), origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularQuadMesh(
    const BaseLib::ISubdivision &div_x,
    const BaseLib::ISubdivision &div_y,
    MathLib::Point3d const& origin,
    std::string const& mesh_name)
{
    std::vector<double> vec_x(div_x());
    std::vector<double> vec_y(div_y());
    std::vector<Node*> nodes(generateRegularNodes(vec_x, vec_y, origin));
    const unsigned n_x_nodes (vec_x.size());

    //elements
    std::vector<Element*> elements;
    const unsigned n_x_cells (vec_x.size()-1);
    const unsigned n_y_cells (vec_y.size()-1);
    elements.reserve(n_x_cells * n_y_cells);

    for (std::size_t j = 0; j < n_y_cells; j++)
    {
        const std::size_t offset_y1 = j * n_x_nodes;
        const std::size_t offset_y2 = (j + 1) * n_x_nodes;
        for (std::size_t k = 0; k < n_x_cells; k++)
        {
            std::array<Node*, 4> element_nodes;
            element_nodes[0] = nodes[offset_y1 + k];
            element_nodes[1] = nodes[offset_y1 + k + 1];
            element_nodes[2] = nodes[offset_y2 + k + 1];
            element_nodes[3] = nodes[offset_y2 + k];
            elements.push_back (new Quad(element_nodes));
        }
    }

    return new Mesh(mesh_name, nodes, elements);
}

Mesh* MeshGenerator::generateRegularHexMesh(
    const double length,
    const std::size_t subdivision,
    const MathLib::Point3d& origin,
    std::string const& mesh_name)
{
    return MeshGenerator::generateRegularHexMesh(subdivision, subdivision,
        subdivision, length/subdivision, origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularHexMesh(
    const double x_length,
    const double y_length,
    const double z_length,
    const std::size_t x_subdivision,
    const std::size_t y_subdivision,
    const std::size_t z_subdivision,
    const MathLib::Point3d& origin,
    std::string const& mesh_name)
{
    return MeshGenerator::generateRegularHexMesh(x_subdivision, y_subdivision, z_subdivision,
        x_length/x_subdivision, y_length/y_subdivision, z_length/z_subdivision, origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularHexMesh(
    const unsigned n_x_cells,
    const unsigned n_y_cells,
    const unsigned n_z_cells,
    const double   cell_size,
    MathLib::Point3d const& origin,
    std::string   const& mesh_name)
{
    return MeshGenerator::generateRegularHexMesh(n_x_cells, n_y_cells, n_z_cells,
        cell_size, cell_size, cell_size, origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularHexMesh(
    const unsigned n_x_cells,
    const unsigned n_y_cells,
    const unsigned n_z_cells,
    const double   cell_size_x,
    const double   cell_size_y,
    const double   cell_size_z,
    MathLib::Point3d const& origin,
    std::string   const& mesh_name)
{
    return generateRegularHexMesh(
            BaseLib::UniformSubdivision(n_x_cells*cell_size_x, n_x_cells),
            BaseLib::UniformSubdivision(n_y_cells*cell_size_y, n_y_cells),
            BaseLib::UniformSubdivision(n_z_cells*cell_size_z, n_z_cells),
            origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularHexMesh(
    const BaseLib::ISubdivision &div_x,
    const BaseLib::ISubdivision &div_y,
    const BaseLib::ISubdivision &div_z,
    MathLib::Point3d const& origin,
    std::string const& mesh_name)
{
    std::vector<double> vec_x(div_x());
    std::vector<double> vec_y(div_y());
    std::vector<double> vec_z(div_z());
    std::vector<Node*> nodes(generateRegularNodes(vec_x, vec_y, vec_z, origin));

    const unsigned n_x_nodes (vec_x.size());
    const unsigned n_y_nodes (vec_y.size());
    const unsigned n_x_cells (vec_x.size()-1);
    const unsigned n_y_cells (vec_y.size()-1);
    const unsigned n_z_cells (vec_z.size()-1);

    //elements
    std::vector<Element*> elements;
    elements.reserve(n_x_cells * n_y_cells * n_z_cells);

    for (std::size_t i = 0; i < n_z_cells; i++)
    {
        const std::size_t offset_z1 = i * n_x_nodes * n_y_nodes; // bottom
        const std::size_t offset_z2 = (i + 1) * n_x_nodes * n_y_nodes; // top
        for (std::size_t j = 0; j < n_y_cells; j++)
        {
            const std::size_t offset_y1 = j * n_x_nodes;
            const std::size_t offset_y2 = (j + 1) * n_x_nodes;
            for (std::size_t k = 0; k < n_x_cells; k++)
            {
                std::array<Node*, 8> element_nodes;
                // bottom
                element_nodes[0] = nodes[offset_z1 + offset_y1 + k];
                element_nodes[1] = nodes[offset_z1 + offset_y1 + k + 1];
                element_nodes[2] = nodes[offset_z1 + offset_y2 + k + 1];
                element_nodes[3] = nodes[offset_z1 + offset_y2 + k];
                // top
                element_nodes[4] = nodes[offset_z2 + offset_y1 + k];
                element_nodes[5] = nodes[offset_z2 + offset_y1 + k + 1];
                element_nodes[6] = nodes[offset_z2 + offset_y2 + k + 1];
                element_nodes[7] = nodes[offset_z2 + offset_y2 + k];
                elements.push_back (new Hex(element_nodes));
            }
        }
    }

    return new Mesh(mesh_name, nodes, elements);
}

Mesh* MeshGenerator::generateRegularTriMesh(
    const double length,
    const std::size_t subdivision,
    const MathLib::Point3d& origin,
    std::string const& mesh_name)
{
    return generateRegularTriMesh(subdivision, subdivision, length/subdivision, origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularTriMesh(
    const double x_length,
    const double y_length,
    const std::size_t x_subdivision,
    const std::size_t y_subdivision,
    const MathLib::Point3d& origin,
    std::string const& mesh_name)
{
    return generateRegularTriMesh(x_subdivision, y_subdivision, x_length/x_subdivision, y_length/y_subdivision, origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularTriMesh(
    const unsigned n_x_cells,
    const unsigned n_y_cells,
    const double cell_size,
    MathLib::Point3d const& origin,
    std::string const& mesh_name)
{
    return generateRegularTriMesh(n_x_cells, n_y_cells, cell_size, cell_size, origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularTriMesh(
    const unsigned n_x_cells,
    const unsigned n_y_cells,
    const double   cell_size_x,
    const double   cell_size_y,
    MathLib::Point3d const& origin,
    std::string   const& mesh_name)
{
    return generateRegularTriMesh(BaseLib::UniformSubdivision(n_x_cells*cell_size_x, n_x_cells),
        BaseLib::UniformSubdivision(n_y_cells*cell_size_y, n_y_cells), origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularTriMesh(
    const BaseLib::ISubdivision &div_x,
    const BaseLib::ISubdivision &div_y,
    MathLib::Point3d const& origin,
    std::string const& mesh_name)
{
    std::vector<double> vec_x(div_x());
    std::vector<double> vec_y(div_y());
    std::vector<Node*> nodes(generateRegularNodes(vec_x, vec_y, origin));
    const unsigned n_x_nodes (vec_x.size());
    const unsigned n_x_cells (vec_x.size()-1);
    const unsigned n_y_cells (vec_y.size()-1);

    //elements
    std::vector<Element*> elements;
    elements.reserve(n_x_cells * n_y_cells * 2);

    for (std::size_t j = 0; j < n_y_cells; j++)
    {
        const std::size_t offset_y1 = j * n_x_nodes;
        const std::size_t offset_y2 = (j + 1) * n_x_nodes;
        for (std::size_t k = 0; k < n_x_cells; k++)
        {
            std::array<Node*, 3> element1_nodes;
            element1_nodes[0] = nodes[offset_y1 + k];
            element1_nodes[1] = nodes[offset_y2 + k + 1];
            element1_nodes[2] = nodes[offset_y2 + k];
            elements.push_back (new Tri(element1_nodes));
            std::array<Node*, 3> element2_nodes;
            element2_nodes[0] = nodes[offset_y1 + k];
            element2_nodes[1] = nodes[offset_y1 + k + 1];
            element2_nodes[2] = nodes[offset_y2 + k + 1];
            elements.push_back (new Tri(element2_nodes));
        }
    }

    return new Mesh(mesh_name, nodes, elements);
}

Mesh* MeshGenerator::generateRegularPrismMesh(
    const double x_length,
    const double y_length,
    const double z_length,
    const std::size_t x_subdivision,
    const std::size_t y_subdivision,
    const std::size_t z_subdivision,
    const MathLib::Point3d& origin,
    std::string const& mesh_name)
{
    return generateRegularPrismMesh(x_subdivision, y_subdivision, z_subdivision,
        x_length/x_subdivision, y_length/y_subdivision, z_length/z_subdivision,
        origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularPrismMesh(
    const unsigned n_x_cells,
    const unsigned n_y_cells,
    const unsigned n_z_cells,
    const double cell_size,
    MathLib::Point3d const& origin,
    std::string const& mesh_name)
{
    return generateRegularPrismMesh(n_x_cells, n_y_cells, n_z_cells,
        cell_size, cell_size, cell_size, origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularPrismMesh(
    const unsigned n_x_cells,
    const unsigned n_y_cells,
    const unsigned n_z_cells,
    const double   cell_size_x,
    const double   cell_size_y,
    const double   cell_size_z,
    MathLib::Point3d const& origin,
    std::string   const& mesh_name)
{
    std::unique_ptr<MeshLib::Mesh> mesh (
        generateRegularTriMesh(n_x_cells, n_y_cells, cell_size_x, cell_size_y, origin, mesh_name));
    std::size_t const n_tris (mesh->getNumberOfElements());
    for (std::size_t i=0; i<n_z_cells; ++i)
        mesh.reset(MeshLib::addTopLayerToMesh(*mesh, cell_size_z, mesh_name));
    std::vector<std::size_t> elem_ids (n_tris);
    std::iota(elem_ids.begin(), elem_ids.end(), 0);
    return MeshLib::removeElements(*mesh, elem_ids, mesh_name);
}

Mesh* MeshGenerator::generateRegularTetMesh(
    const double x_length,
    const double y_length,
    const double z_length,
    const std::size_t x_subdivision,
    const std::size_t y_subdivision,
    const std::size_t z_subdivision,
    const MathLib::Point3d& origin,
    std::string const& mesh_name)
{
    return generateRegularTetMesh(x_subdivision, y_subdivision, z_subdivision,
        x_length/x_subdivision, y_length/y_subdivision, z_length/z_subdivision,
        origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularTetMesh(
    const unsigned n_x_cells,
    const unsigned n_y_cells,
    const unsigned n_z_cells,
    const double   cell_size_x,
    const double   cell_size_y,
    const double   cell_size_z,
    MathLib::Point3d const& origin,
    std::string   const& mesh_name)
{
    return generateRegularTetMesh(
            BaseLib::UniformSubdivision(n_x_cells*cell_size_x, n_x_cells),
            BaseLib::UniformSubdivision(n_y_cells*cell_size_y, n_y_cells),
            BaseLib::UniformSubdivision(n_z_cells*cell_size_z, n_z_cells),
            origin, mesh_name);
}

Mesh* MeshGenerator::generateRegularTetMesh(
    const BaseLib::ISubdivision &div_x,
    const BaseLib::ISubdivision &div_y,
    const BaseLib::ISubdivision &div_z,
    MathLib::Point3d const& origin,
    std::string const& mesh_name)
{
    std::vector<double> vec_x(div_x());
    std::vector<double> vec_y(div_y());
    std::vector<double> vec_z(div_z());
    std::vector<Node*> nodes(generateRegularNodes(vec_x, vec_y, vec_z, origin));

    const unsigned n_x_nodes (vec_x.size());
    const unsigned n_y_nodes (vec_y.size());
    const unsigned n_x_cells (vec_x.size()-1);
    const unsigned n_y_cells (vec_y.size()-1);
    const unsigned n_z_cells (vec_z.size()-1);

    //elements
    std::vector<Element*> elements;
    elements.reserve(n_x_cells * n_y_cells * n_z_cells * 5);

    for (std::size_t i = 0; i < n_z_cells; i++)
    {
        const std::size_t offset_z1 = i * n_x_nodes * n_y_nodes; // bottom
        const std::size_t offset_z2 = (i + 1) * n_x_nodes * n_y_nodes; // top
        for (std::size_t j = 0; j < n_y_cells; j++)
        {
            const std::size_t offset_y1 = j * n_x_nodes;
            const std::size_t offset_y2 = (j + 1) * n_x_nodes;
            for (std::size_t k = 0; k < n_x_cells; k++)
            {
                // tet 1
                {
                    std::array<Node*, 4> element_nodes;
                    // bottom
                    element_nodes[0] = nodes[offset_z1 + offset_y1 + k];
                    element_nodes[1] = nodes[offset_z1 + offset_y1 + k + 1];
                    element_nodes[2] = nodes[offset_z1 + offset_y2 + k + 1];
                    // top
                    element_nodes[3] = nodes[offset_z2 + offset_y1 + k + 1];
                    elements.push_back (new Tet(element_nodes));
                }
                // tet 2
                {
                    std::array<Node*, 4> element_nodes;
                    // bottom
                    element_nodes[0] = nodes[offset_z1 + offset_y1 + k];
                    element_nodes[1] = nodes[offset_z1 + offset_y2 + k + 1];
                    element_nodes[2] = nodes[offset_z1 + offset_y2 + k];
                    // top
                    element_nodes[3] = nodes[offset_z2 + offset_y2 + k];
                    elements.push_back (new Tet(element_nodes));
                }
                // tet 3
                {
                    std::array<Node*, 4> element_nodes;
                    // bottom
                    element_nodes[0] = nodes[offset_z1 + offset_y1 + k];
                    // top
                    element_nodes[1] = nodes[offset_z2 + offset_y1 + k];
                    element_nodes[2] = nodes[offset_z2 + offset_y1 + k + 1];
                    element_nodes[3] = nodes[offset_z2 + offset_y2 + k];
                    elements.push_back (new Tet(element_nodes));
                }
                // tet 4
                {
                    std::array<Node*, 4> element_nodes;
                    // bottom
                    element_nodes[0] = nodes[offset_z1 + offset_y2 + k + 1];
                    // top
                    element_nodes[1] = nodes[offset_z2 + offset_y1 + k + 1];
                    element_nodes[2] = nodes[offset_z2 + offset_y2 + k + 1];
                    element_nodes[3] = nodes[offset_z2 + offset_y2 + k];
                    elements.push_back (new Tet(element_nodes));
                }
                // tet 5
                {
                    std::array<Node*, 4> element_nodes;
                    // bottom
                    element_nodes[0] = nodes[offset_z1 + offset_y1 + k];
                    element_nodes[1] = nodes[offset_z1 + offset_y2 + k + 1];
                    // top
                    element_nodes[2] = nodes[offset_z2 + offset_y2 + k];
                    element_nodes[3] = nodes[offset_z2 + offset_y1 + k + 1];
                    elements.push_back (new Tet(element_nodes));
                }
            }
        }
    }

    return new Mesh(mesh_name, nodes, elements);
}

MeshLib::Mesh*
MeshGenerator::createSurfaceMesh(std::string const& mesh_name,
    MathLib::Point3d const& ll, MathLib::Point3d const& ur,
    std::array<std::size_t, 2> const& n_steps,
    const std::function<double(double,double)>& f)
{
    std::array<double, 2> step_size{{
        (ur[0]-ll[0])/(n_steps[0]-1), (ur[1]-ll[1])/(n_steps[1]-1)}};

    std::vector<MeshLib::Node*> nodes;
    for (std::size_t j(0); j<n_steps[1]; ++j) {
        for (std::size_t i(0); i<n_steps[0]; ++i) {
            std::size_t const id(i+j*n_steps[1]);
            std::array<double, 3> coords;
            coords[0] = ll[0]+i*step_size[0];
            coords[1] = ll[1]+j*step_size[1];
            coords[2] = f(coords[0],coords[1]);
            nodes.push_back(new MeshLib::Node(coords, id));
        }
    }

    std::vector<MeshLib::Element*> sfc_eles;
    for (std::size_t j(0); j<n_steps[1]-1; ++j) {
        for (std::size_t i(0); i<n_steps[0]-1; ++i) {
            std::size_t id_ll(i+j*n_steps[0]);
            std::size_t id_lr(i+1+j*n_steps[0]);
            std::size_t id_ul(i+(j+1)*n_steps[0]);
            std::size_t id_ur(i+1+(j+1)*n_steps[0]);
            sfc_eles.push_back(new MeshLib::Tri(std::array<MeshLib::Node*,3>
                {{nodes[id_ll], nodes[id_lr], nodes[id_ur]}}));
            sfc_eles.push_back(new MeshLib::Tri(std::array<MeshLib::Node*,3>
                {{nodes[id_ll], nodes[id_ur], nodes[id_ul]}}));
        }
    }

    return new MeshLib::Mesh(mesh_name, nodes, sfc_eles);
}

}
