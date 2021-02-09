/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

// ThirdParty
#include <tclap/CmdLine.h>

#include "GeoLib/AABB.h"
#include "InfoLib/GitInfo.h"
#include "MathLib/Point3d.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

static std::string mat_name = "MaterialIDs";

// tests if point p is located outside of given AABB
bool testPointOutsideAABB(MathLib::Point3d const& p,
                          MathLib::Point3d const& min_pnt,
                          MathLib::Point3d const& max_pnt)
{
    for (std::size_t i = 0; i < 3; ++i)
    {
        if (p[i] < min_pnt[i] || p[i] > max_pnt[i])
        {
            return true;
        }
    }
    return false;
}

// tests if a plane and an AABB are intersecting
// (based on Christer Ericson "Real Time Collision Detection" 5.2.3)
bool testAABBIntersectingPlane(Eigen::Vector3d const& aabb_centre,
                               Eigen::Vector3d const& aabb_extent,
                               Eigen::Vector3d const& plane_normal,
                               double const pd)
{

    double const r = aabb_extent.dot(plane_normal.cwiseAbs());
    double const s = plane_normal.dot(aabb_centre) - pd;
    return std::abs(s) > r;
}

// separating axes test for triangle vs AABB projection
// (based on Christer Ericson "Real Time Collision Detection" 5.2.9)
bool paxx(Eigen::Vector3d const& v0, Eigen::Vector3d const& v1,
          Eigen::Vector3d const& v2, Eigen::Vector3d const& a,
          Eigen::Vector3d const& e)
{
    double const p0 = v0.dot(a);
    double const p1 = v1.dot(a);
    double const p2 = v2.dot(a);
    double const r = e.x() * std::abs(a.x()) + e.y() * std::abs(a.y()) +
                     e.z() * std::abs(a.z());
    return std::max(-std::max({p0, p1, p2}), std::min({p0, p1, p2})) <= r;
}

// tests if a triangle and an AABB are intersecting
// (based on Christer Ericson "Real Time Collision Detection" 5.2.9)
bool testTriangleIntersectingAABB(MeshLib::Node const& n0,
                                  MeshLib::Node const& n1,
                                  MeshLib::Node const& n2,
                                  Eigen::Vector3d const& c,
                                  Eigen::Vector3d const& e)
{
    // Translate triangle as conceptually moving AABB to origin
    Eigen::Vector3d const v0(Eigen::Vector3d(n0.getCoords()) - c);
    Eigen::Vector3d const v1(Eigen::Vector3d(n1.getCoords()) - c);
    Eigen::Vector3d const v2(Eigen::Vector3d(n2.getCoords()) - c);

    // Test the three axes corresponding to the face normals of AABB b
    if (std::max({v0.x(), v1.x(), v2.x()}) < -e.x() ||
        std::min({v0.x(), v1.x(), v2.x()}) > e.x())
    {
        return false;
    }
    if (std::max({v0.y(), v1.y(), v2.y()}) < -e.y() ||
        std::min({v0.y(), v1.y(), v2.y()}) > e.y())
    {
        return false;
    }
    if (std::max({v0.z(), v1.z(), v2.z()}) < -e.z() ||
        std::min({v0.z(), v1.z(), v2.z()}) > e.z())
    {
        return false;
    }

    // Compute edge vectors for triangle
    Eigen::Vector3d const tri_edge0(v1 - v0);
    Eigen::Vector3d const tri_edge1(v2 - v1);
    Eigen::Vector3d const tri_edge2(v0 - v2);

    // separating axes
    std::array<Eigen::Vector3d, 9> axx;
    axx[0] = Eigen::Vector3d({0, -tri_edge0.z(), tri_edge0.y()});
    axx[1] = Eigen::Vector3d({0, -tri_edge1.z(), tri_edge1.y()});
    axx[2] = Eigen::Vector3d({0, -tri_edge2.z(), tri_edge2.y()});
    axx[3] = Eigen::Vector3d({tri_edge0.z(), 0, -tri_edge0.x()});
    axx[4] = Eigen::Vector3d({tri_edge1.z(), 0, -tri_edge1.x()});
    axx[5] = Eigen::Vector3d({tri_edge2.z(), 0, -tri_edge2.x()});
    axx[6] = Eigen::Vector3d({-tri_edge0.y(), tri_edge0.x(), 0});
    axx[7] = Eigen::Vector3d({-tri_edge1.y(), tri_edge1.x(), 0});
    axx[8] = Eigen::Vector3d({-tri_edge2.y(), tri_edge2.x(), 0});

    // testing separating axes against triangle
    for (auto const& a : axx)
    {
        if (!paxx(v0, v1, v2, a, e))
            return false;
    }

    // Test separating axis corresponding to triangle face normal
    Eigen::Vector3d const plane_normal(tri_edge0.cross(tri_edge1));
    double const pd = plane_normal.dot(v0);
    return testAABBIntersectingPlane(c, e, plane_normal, pd);
}

void markFaults(MeshLib::Mesh& mesh, MeshLib::Mesh const& fault,
                int const fault_id, std::array<double, 3> half_cell_size)
{
    auto const& elems = mesh.getElements();
    std::size_t const n_elems = mesh.getNumberOfElements();
    auto mat_ids = mesh.getProperties().getPropertyVector<int>(mat_name);
    auto const& fnodes = fault.getNodes();
    auto const& felems = fault.getElements();
    std::size_t const n_felems = fault.getNumberOfElements();
    GeoLib::AABB fault_aabb(fnodes.cbegin(), fnodes.cend());
    auto min_pnt = fault_aabb.getMinPoint();
    auto max_pnt = fault_aabb.getMaxPoint();

    // get bounding box of fault + voxel extent
    for (std::size_t i = 0; i < 3; ++i)
    {
        min_pnt[i] -= half_cell_size[i];
        max_pnt[i] += half_cell_size[i];
    }

    // test each voxel grid element vs each fault triangle
    Eigen::Vector3d const extent(
        {half_cell_size[0], half_cell_size[1], half_cell_size[2]});
    for (std::size_t j = 0; j < n_elems; ++j)
    {
        // test if bounding box of fault is intersecting voxel
        auto const& centre_pnt = MeshLib::getCenterOfGravity(*elems[j]);
        if (testPointOutsideAABB(centre_pnt, min_pnt, max_pnt))
        {
            continue;
        }

        // test if voxel is intersecting triangle
        Eigen::Vector3d const c(centre_pnt.getCoords());
        for (auto const* const fault_elem : felems)
        {
            if (fault_elem->getDimension() != 2)
            {
                continue;
            }

            if (testTriangleIntersectingAABB(
                    *fault_elem->getNode(0), *fault_elem->getNode(1),
                    *fault_elem->getNode(2), c, extent))
            {
                (*mat_ids)[j] = fault_id;
                break;
            }

            if (fault_elem->getGeomType() == MeshLib::MeshElemType::QUAD)
            {
                if (testTriangleIntersectingAABB(
                        *fault_elem->getNode(0), *fault_elem->getNode(2),
                        *fault_elem->getNode(3), c, extent))
                {
                    (*mat_ids)[j] = fault_id;
                    break;
                }
            }
        }
    }
}

// test if input mesh is voxel grid
bool isVoxelGrid(MeshLib::Mesh const& mesh)
{
    auto const& elements = mesh.getElements();
    if (std::any_of(elements.cbegin(), elements.cend(), [&](auto const& e) {
            return (e->getGeomType() != MeshLib::MeshElemType::HEXAHEDRON);
        }))
    {
        ERR("Input mesh needs to be voxel grid (i.e. equally sized axis "
            "aligned hexahedra).");
        return false;
    }

    for (auto const& e : elements)
    {
        auto const n = e->getNodes();
        if ((*n[0])[2] != (*n[1])[2] || (*n[1])[2] != (*n[2])[2] ||
            (*n[4])[2] != (*n[5])[2] || (*n[5])[2] != (*n[6])[2] ||
            (*n[1])[0] != (*n[2])[0] || (*n[2])[0] != (*n[5])[0] ||
            (*n[0])[0] != (*n[3])[0] || (*n[3])[0] != (*n[7])[0])
        {
            ERR("Input mesh needs to be voxel grid (i.e. equally sized axis "
                "aligned hexahedra).");
            return false;
        }
    }
    return true;
}

int main(int argc, char* argv[])
{
    constexpr int mat_not_set = std::numeric_limits<int>::max();

    TCLAP::CmdLine cmd(
        "Marks all elements in a voxel grid (i.e. a structured hex grid, for "
        "instance created with Layers2Grid or Vtu2Grid) that are intersected "
        "by a triangulated 2D mesh representing a fault or some other "
        "significant structure. The material group for those intersected "
        "elements can be explicitely specified, otherwise the largest existing "
        "MaterialID will be increased by one.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2021, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<int> id_arg(
        "m", "material", "material id for cells intersected by fault", false,
        mat_not_set, "non-negative integer");
    cmd.add(id_arg);

    TCLAP::ValueArg<std::string> output_arg(
        "o", "output", "name of output mesh (*.vtu)", true, "", "string");
    cmd.add(output_arg);

    TCLAP::ValueArg<std::string> fault_arg(
        "f", "fault", "name of mesh representing fault (*.vtu)", true, "",
        "string");
    cmd.add(fault_arg);

    TCLAP::ValueArg<std::string> input_arg(
        "i", "input",
        "name of the input file list containing the paths the all input layers "
        "in correct order from top to bottom",
        true, "", "string");
    cmd.add(input_arg);
    cmd.parse(argc, argv);

    std::string const input_name = input_arg.getValue();
    std::string const fault_name = fault_arg.getValue();
    std::string const output_name = output_arg.getValue();

    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::readMeshFromFile(input_name));
    if (!isVoxelGrid(*mesh))
    {
        return EXIT_FAILURE;
    }
    auto const& mat_ids = MeshLib::materialIDs(*mesh);
    if (!mat_ids)
    {
        ERR("Input mesh has no material IDs");
        return EXIT_FAILURE;
    }

    std::unique_ptr<MeshLib::Mesh> fault(
        MeshLib::IO::readMeshFromFile(fault_name));
    if (fault->getDimension() != 2)
    {
        ERR("Fault needs to be a 2D mesh.");
        return EXIT_FAILURE;
    }

    int fault_id = id_arg.getValue();
    if (!id_arg.isSet())
    {
        auto it = std::max_element(mat_ids->cbegin(), mat_ids->cend());
        fault_id = *it + 1;
    }

    std::array<double, 3> half_cell_size;
    {
        auto const n = *mesh->getElement(0)->getNode(0);
        auto const c = MeshLib::getCenterOfGravity(*mesh->getElement(0));
        half_cell_size[0] = std::abs(c[0] - n[0]);
        half_cell_size[1] = std::abs(c[1] - n[1]);
        half_cell_size[2] = std::abs(c[2] - n[2]);
    };

    markFaults(*mesh, *fault, fault_id, half_cell_size);

    MeshLib::IO::VtuInterface vtu(mesh.get());
    vtu.writeToFile(output_name);
    return EXIT_SUCCESS;
}
