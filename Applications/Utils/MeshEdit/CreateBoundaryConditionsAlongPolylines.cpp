/*
 * \date 2014-09-30
 * \brief Create BoundaryConditions from a polylines.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <map>
#include <string>
#include <vector>
#include <fstream>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/FileTools.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "GeoLib/IO/readGeometryFromFile.h"
#include "GeoLib/IO/writeGeometryToFile.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/MeshSurfaceExtraction.h"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"

void convertMeshNodesToGeometry(std::vector<MeshLib::Node*> const& nodes,
    std::vector<std::size_t> const& node_ids,
    std::string & geo_name,
    GeoLib::GEOObjects & geometry_sets)
{
    // copy data
    auto pnts = std::unique_ptr<std::vector<GeoLib::Point*>>(
        new std::vector<GeoLib::Point*>);
    std::unique_ptr<std::map<std::string, std::size_t>> pnt_names{
        new std::map<std::string, std::size_t>};
    std::size_t cnt(0);
    for (std::size_t id: node_ids) {
        pnts->push_back(new GeoLib::Point(*(nodes[id]), cnt));
        pnt_names->insert(std::pair<std::string, std::size_t>(
            geo_name+"-PNT-"+std::to_string(cnt), cnt));
        cnt++;
    }

    // create data structures for geometry
    geometry_sets.addPointVec(std::move(pnts), geo_name, std::move(pnt_names));
}

void writeGroundwaterFlowPointBC(std::ostream& bc_out,
                                 std::string const& pnt_name, double head_value)
{
    bc_out << "#BOUNDARY_CONDITION\n";
    bc_out << "  $PCS_TYPE\n";
    bc_out << "    GROUNDWATER_FLOW\n";
    bc_out << "  $PRIMARY_VARIABLE\n";
    bc_out << "    HEAD\n";
    bc_out << "  $GEO_TYPE\n";
    bc_out << "    POINT " << pnt_name << "\n";
    bc_out << "  $DIS_TYPE\n";
    bc_out << "    CONSTANT " << head_value << "\n";
}

void writeLiquidFlowPointBC(std::ostream & bc_out, std::string const& pnt_name)
{
    bc_out << "#BOUNDARY_CONDITION\n";
    bc_out << "  $PCS_TYPE\n";
    bc_out << "    LIQUID_FLOW\n";
    bc_out << "  $PRIMARY_VARIABLE\n";
    bc_out << "    PRESSURE1\n";
    bc_out << "  $GEO_TYPE\n";
    bc_out << "    POINT " << pnt_name << "\n";
    bc_out << "  $DIS_TYPE\n";
    bc_out << "    CONSTANT 0.0\n";
}

// geometry_sets contains the geometric points the boundary conditions will be
// set on, geo_name is the name the geometry can be accessed with, out_fname is
// the base file name the gli and bc as well as the gml file will be written to.
void writeBCsAndGeometry(GeoLib::GEOObjects& geometry_sets,
                         std::string& geo_name, std::string const& out_fname,
                         std::string const& bc_type, bool write_gml)
{
    if (write_gml) {
        INFO("write points to \"%s.gml\".", geo_name.c_str());
        GeoLib::IO::writeGeometryToFile(geo_name, geometry_sets, out_fname+".gml");
    }
    GeoLib::IO::writeGeometryToFile(geo_name, geometry_sets, out_fname+".gli");

    bool liquid_flow(false);
    if (bc_type == "LIQUID_FLOW")
        liquid_flow = true;


    GeoLib::PointVec const* pnt_vec_objs(geometry_sets.getPointVecObj(geo_name));
    std::vector<GeoLib::Point*> const& pnts(*(pnt_vec_objs->getVector()));
    std::ofstream bc_out (out_fname+".bc");
    for (std::size_t k(0); k<pnts.size(); k++) {
        std::string const& pnt_name(pnt_vec_objs->getItemNameByID(k));
        if (!pnt_name.empty()) {
            if (liquid_flow)
                writeLiquidFlowPointBC(bc_out, pnt_name);
            else
                writeGroundwaterFlowPointBC(bc_out, pnt_name, (*pnts[k])[2]);
        }
    }
    bc_out << "#STOP\n";
    bc_out.close();
}

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd(
        "Creates boundary conditions for mesh nodes along polylines."
        "The documentation is available at https://docs.opengeosys.org/docs/tools/model-preparation/create-boundary-conditions-along-a-polyline",
        ' ',
        "0.1");
    TCLAP::ValueArg<bool> gml_arg("", "gml",
        "if switched on write found nodes to file in gml format", false, false, "bool");
    cmd.add(gml_arg);

    TCLAP::ValueArg<std::string> output_base_fname("o", "output-base-file-name",
        "the base name of the file the output (geometry (gli) and boundary"\
        "condition (bc)) will be written to", true,
        "", "file name");
    cmd.add(output_base_fname);

    TCLAP::ValueArg<std::string> bc_type("t", "type",
        "the process type the boundary condition will be written for "\
        "currently LIQUID_FLOW (primary variable PRESSURE1) and "\
        "GROUNDWATER_FLOW (primary variable HEAD, default) are supported", true,
        "",
        "process type as string (LIQUID_FLOW or GROUNDWATER_FLOW (default))");
    cmd.add(bc_type);

    TCLAP::ValueArg<double> search_length_arg("s", "search-length",
        "The size of the search length. The default value is "
        "std::numeric_limits<double>::epsilon()", false,
        std::numeric_limits<double>::epsilon(), "floating point number");
    cmd.add(search_length_arg);

    TCLAP::ValueArg<std::string> geometry_fname("i", "input-geometry",
        "the name of the file containing the input geometry", true,
        "", "file name");
    cmd.add(geometry_fname);

    TCLAP::ValueArg<std::string> mesh_arg("m", "mesh-file",
        "the name of the file containing the mesh", true,
        "", "file name");
    cmd.add(mesh_arg);

    cmd.parse(argc, argv);

    // *** read mesh
    INFO("Reading mesh \"%s\" ... ", mesh_arg.getValue().c_str());
    std::unique_ptr<MeshLib::Mesh> subsurface_mesh(
        MeshLib::IO::readMeshFromFile(mesh_arg.getValue()));
    INFO("done.");
    INFO("Extracting top surface of mesh \"%s\" ... ",
        mesh_arg.getValue().c_str());
    const MathLib::Vector3 dir(0,0,-1);
    double const angle(90);
    std::unique_ptr<MeshLib::Mesh> surface_mesh(
        MeshLib::MeshSurfaceExtraction::getMeshSurface(*subsurface_mesh, dir,
                                                       angle));
    INFO("done.");
    subsurface_mesh.reset(nullptr);

    // *** read geometry
    GeoLib::GEOObjects geometries;
    GeoLib::IO::readGeometryFromFile(geometry_fname.getValue(), geometries);

    std::string geo_name;
    {
        std::vector<std::string> geo_names;
        geometries.getGeometryNames(geo_names);
        geo_name = geo_names[0];
    }

    // *** check if the data is usable
    // *** get vector of polylines
    std::vector<GeoLib::Polyline*> const* plys(geometries.getPolylineVec(geo_name));
    if (!plys) {
        ERR("Could not get vector of polylines out of geometry \"%s\".",
            geo_name.c_str());
        return EXIT_FAILURE;
    }

    MeshGeoToolsLib::SearchLength search_length_strategy;
    if (search_length_arg.isSet()) {
        search_length_strategy =
            MeshGeoToolsLib::SearchLength(search_length_arg.getValue());
    }

    GeoLib::GEOObjects geometry_sets;
    MeshGeoToolsLib::MeshNodeSearcher mesh_searcher(
        *surface_mesh, std::move(search_length_strategy),
        MeshGeoToolsLib::SearchAllNodes::Yes);
    for(std::size_t k(0); k<plys->size(); k++) {
        std::vector<std::size_t> ids
            (mesh_searcher.getMeshNodeIDsAlongPolyline(*((*plys)[k])));
        if (ids.empty())
            continue;
        std::string geo_name("Polyline-"+std::to_string(k));
        convertMeshNodesToGeometry(surface_mesh->getNodes(), ids, geo_name,
            geometry_sets);
    }

    // merge all together
    std::vector<std::string> geo_names;
    geometry_sets.getGeometryNames(geo_names);
    if (geo_names.empty()) {
        ERR("Did not find mesh nodes along polylines.");
        return EXIT_FAILURE;
    }

    std::string merge_name("AllMeshNodesAlongPolylines");
    if (geometry_sets.mergeGeometries(geo_names, merge_name) == 2)
        merge_name = geo_names[0];

    GeoLib::PointVec const* pnt_vec(geometry_sets.getPointVecObj(merge_name));
    std::vector<GeoLib::Point*> const* merged_pnts(pnt_vec->getVector());

    std::vector<GeoLib::Point> pnts_with_id;
    const std::size_t n_merged_pnts(merged_pnts->size());
    for(std::size_t k(0); k<n_merged_pnts; ++k) {
        pnts_with_id.emplace_back(*((*merged_pnts)[k]), k);
    }

    std::sort(pnts_with_id.begin(), pnts_with_id.end(),
        [](GeoLib::Point const& p0, GeoLib::Point const& p1)
            { return p0 < p1; }
    );

    double const eps (std::numeric_limits<double>::epsilon());
    auto surface_pnts = std::unique_ptr<std::vector<GeoLib::Point*>>(
        new std::vector<GeoLib::Point*>);
    std::unique_ptr<std::map<std::string, std::size_t>> name_id_map{
        new std::map<std::string, std::size_t>};

    // insert first point
    surface_pnts->push_back(
        new GeoLib::Point(pnts_with_id[0], surface_pnts->size()));
    std::string element_name;
    pnt_vec->getNameOfElementByID(0, element_name);
    name_id_map->insert(
        std::pair<std::string, std::size_t>(element_name,0)
    );
    for (std::size_t k(1); k < n_merged_pnts; ++k) {
        const GeoLib::Point& p0 (pnts_with_id[k-1]);
        const GeoLib::Point& p1 (pnts_with_id[k]);
        if (std::abs (p0[0] - p1[0]) > eps || std::abs (p0[1] - p1[1]) > eps) {
            surface_pnts->push_back(new GeoLib::Point(pnts_with_id[k],
                surface_pnts->size()));
            std::string element_name;
            pnt_vec->getNameOfElementByID(k, element_name);
            name_id_map->insert(
                std::pair<std::string, std::size_t>(element_name,
                surface_pnts->size()-1)
            );
        }
    }

    std::string surface_name(BaseLib::dropFileExtension(mesh_arg.getValue())+"-MeshNodesAlongPolylines");
    geometry_sets.addPointVec(std::move(surface_pnts), surface_name,
                              std::move(name_id_map), 1e-6);

    // write the BCs and the merged geometry set to file
    std::string const base_fname(
        BaseLib::dropFileExtension(output_base_fname.getValue()));
    writeBCsAndGeometry(geometry_sets, surface_name, base_fname,
                        bc_type.getValue(), gml_arg.getValue());
    return EXIT_SUCCESS;
}
