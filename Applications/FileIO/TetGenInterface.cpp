/**
 * \file
 * \author Thomas Fischer
 * \date   2011-09-12
 * \brief  Implementation of the TetGenInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TetGenInterface.h"

#include <cstddef>
#include <string>
#include <fstream>

#include <logog/include/logog.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "GeoLib/Triangle.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/MeshInformation.h"

namespace FileIO
{
TetGenInterface::TetGenInterface() :
    _zero_based_idx (false), _boundary_markers (false)
{
}

bool TetGenInterface::readTetGenGeometry (std::string const& geo_fname,
                                          GeoLib::GEOObjects &geo_objects)
{
    std::ifstream poly_stream (geo_fname.c_str());

    if (!poly_stream)
    {
        ERR ("TetGenInterface::readTetGenGeometry() failed to open %s", geo_fname.c_str());
        return false;
    }
    std::string ext (BaseLib::getFileExtension(geo_fname));
    if (ext != "smesh")
    {
        ERR ("TetGenInterface::readTetGenGeometry() - unknown file type (only *.smesh is supported).");
        return false;
    }

    std::vector<MeshLib::Node*> nodes;
    if (!readNodesFromStream (poly_stream, nodes))
    {
        // remove nodes read until now
        for (auto & node : nodes)
            delete node;
        return false;
    }
    const std::size_t nNodes (nodes.size());
    auto points = std::make_unique<std::vector<GeoLib::Point*>>();
    points->reserve(nNodes);
    for (std::size_t k(0); k<nNodes; ++k)
    {
        points->push_back(new GeoLib::Point(*(nodes[k]), nodes[k]->getID()));
        delete nodes[k];
    }
    std::string geo_name (BaseLib::extractBaseNameWithoutExtension(geo_fname));
    geo_objects.addPointVec(std::move(points), geo_name);
    const std::vector<std::size_t> &id_map (geo_objects.getPointVecObj(geo_name)->getIDMap());

    auto surfaces = std::make_unique<std::vector<GeoLib::Surface*>>();
    if (!parseSmeshFacets(poly_stream, *surfaces, *geo_objects.getPointVec(geo_name), id_map))
    {
        // remove surfaces read until now but keep the points
        for (std::size_t k=0; k<surfaces->size(); k++)
            delete (*surfaces)[k];
    }
    geo_objects.addSurfaceVec(std::move(surfaces), geo_name);

    return true;
}

std::size_t TetGenInterface::getNFacets(std::ifstream &input)
{
    std::string line;
    while (!input.fail())
    {
        getline (input, line);
        if (input.fail())
        {
            ERR("TetGenInterface::getNFacets(): Error reading number of facets.");
            return 0;
        }

        BaseLib::simplify(line);
        if (line.empty() || line.compare(0,1,"#") == 0)
            continue;

        const std::list<std::string> fields = BaseLib::splitString(line, ' ');
        auto it = fields.begin();
        const auto nFacets(BaseLib::str2number<std::size_t>(*it));
        if (fields.size() > 1)
            _boundary_markers = BaseLib::str2number<std::size_t>(*(++it)) != 0;
        return nFacets;
    }
    return 0;
}

bool TetGenInterface::parseSmeshFacets(std::ifstream &input,
                                       std::vector<GeoLib::Surface*> &surfaces,
                                       const std::vector<GeoLib::Point*> &points,
                                       const std::vector<std::size_t> &pnt_id_map)
{
    const std::size_t nFacets (this->getNFacets(input));
    std::string line;
    surfaces.reserve(nFacets);
    std::list<std::string>::const_iterator it;

    const unsigned offset = (_zero_based_idx) ? 0 : 1;
    std::vector<std::size_t> idx_map;

    std::size_t k(0);
    while (k<nFacets && !input.fail())
    {
        getline (input, line);
        if (input.fail())
        {
            ERR("TetGenInterface::parseFacets(): Error reading facet %d.", k);
            return false;
        }

        BaseLib::simplify(line);
        if (line.empty() || line.compare(0,1,"#") == 0)
        {
            continue;
        }

        // read facets
        const std::list<std::string> point_fields = BaseLib::splitString(line, ' ');
        it = point_fields.begin();
        const auto nPoints = BaseLib::str2number<std::size_t>(*it);
        if (nPoints != 3)
        {
            ERR ("Smesh-files are currently only supported for triangle meshes.");
            return false;
        }
        std::vector<std::size_t> point_ids;
        const std::size_t point_field_size = (_boundary_markers) ? nPoints+1 : nPoints;
        if (point_fields.size() > point_field_size)
        {
            for (std::size_t j(0); j<nPoints; ++j)
                point_ids.push_back(pnt_id_map[BaseLib::str2number<std::size_t>(*(++it))-offset]);

            const std::size_t sfc_marker = (_boundary_markers) ? BaseLib::str2number<std::size_t>(*(++it)) : 0;
            const std::size_t idx = std::find(idx_map.begin(), idx_map.end(), sfc_marker) - idx_map.begin();
            if (idx >= surfaces.size())
            {
                idx_map.push_back(sfc_marker);
                surfaces.push_back(new GeoLib::Surface(points));
            }
            surfaces[idx]->addTriangle(point_ids[0], point_ids[1], point_ids[2]);
        }
        else
        {
            ERR("TetGenInterface::parseFacets(): Error reading points for facet %d.", k);
            return false;
        }
        ++k;
    }
    // here the poly-file potentially defines a number of points to mark holes within the volumes defined by the facets, these are ignored for now
    // here the poly-file potentially defines a number of region attributes, these are ignored for now

    std::size_t nTotalTriangles (0);
    for (auto & surface : surfaces)
        nTotalTriangles += surface->getNumberOfTriangles();
    if (nTotalTriangles == nFacets)
        return true;

    ERR ("TetGenInterface::parseFacets(): Number of expected total triangles (%d) does not match number of found triangles (%d).", surfaces.size(), nTotalTriangles);
    return false;
}

MeshLib::Mesh* TetGenInterface::readTetGenMesh (std::string const& nodes_fname,
                                                std::string const& ele_fname)
{
    std::ifstream ins_nodes (nodes_fname.c_str());
    std::ifstream ins_ele (ele_fname.c_str());

    if (!ins_nodes || !ins_ele)
    {
        if (!ins_nodes)
            ERR ("TetGenInterface::readTetGenMesh failed to open %s", nodes_fname.c_str());
        if (!ins_ele)
            ERR ("TetGenInterface::readTetGenMesh failed to open %s", ele_fname.c_str());
        return nullptr;
    }

    std::vector<MeshLib::Node*> nodes;
    if (!readNodesFromStream (ins_nodes, nodes)) {
        // remove nodes read until now
        for (auto & node : nodes) {
            delete node;
        }
        return nullptr;
    }

    std::vector<MeshLib::Element*> elements;
    std::vector<int> materials;
    if (!readElementsFromStream (ins_ele, elements, materials, nodes)) {
        // remove elements read until now
        for (auto & element : elements) {
            delete element;
        }
        // remove nodes
        for (auto & node : nodes) {
            delete node;
        }
        return nullptr;
    }

    MeshLib::Properties properties;
    // Transmit material values if there is any material value != 0
    if (std::any_of(materials.cbegin(), materials.cend(), [](int m){ return m != 0; }))
    {
        auto* const mat_props = properties.createNewPropertyVector<int>(
            "MaterialIDs", MeshLib::MeshItemType::Cell);
        mat_props->reserve(elements.size());
        std::copy(materials.cbegin(),
                  materials.cend(),
                  std::back_inserter(*mat_props));
    }

    const std::string mesh_name (BaseLib::extractBaseNameWithoutExtension(nodes_fname));
    return new MeshLib::Mesh(mesh_name, nodes, elements, properties);
}

bool TetGenInterface::readNodesFromStream (std::ifstream &ins,
                                           std::vector<MeshLib::Node*> &nodes)
{
    std::string line;
    getline (ins, line);
    std::size_t n_nodes, dim, n_attributes;
    bool boundary_markers;

    while (!ins.fail())
    {
        BaseLib::simplify(line);
        if (line.empty() || line.compare(0,1,"#") == 0)
        {
            // this line is a comment - skip
            getline (ins, line);
            continue;
        }
        // read header line
        bool header_okay = parseNodesFileHeader(line, n_nodes, dim, n_attributes, boundary_markers);
        if (!header_okay)
            return false;
        if (!parseNodes(ins, nodes, n_nodes, dim))
            return false;
        return true;
    }
    return false;
}

bool TetGenInterface::parseNodesFileHeader(std::string &line,
                                           std::size_t &n_nodes,
                                           std::size_t &dim,
                                           std::size_t &n_attributes,
                                           bool &boundary_markers) const
{
    std::list<std::string> pnt_header = BaseLib::splitString(line, ' ');
    if (pnt_header.empty())
    {
        ERR("TetGenInterface::parseNodesFileHeader(): could not read number of nodes specified in header.");
        return false;
    }
    auto it = pnt_header.begin();
    n_nodes = BaseLib::str2number<std::size_t> (*it);
    dim = (pnt_header.size()==1) ? 3 : BaseLib::str2number<std::size_t> (*(++it));

    if (pnt_header.size()<4)
    {
        n_attributes = 0;
        boundary_markers = false;
        return true;
    }

    n_attributes = BaseLib::str2number<std::size_t> (*(++it));
    boundary_markers = *(++it) == "1";

    return true;
}

bool TetGenInterface::parseNodes(std::ifstream &ins,
                                 std::vector<MeshLib::Node*> &nodes,
                                 std::size_t n_nodes,
                                 std::size_t dim)
{
    std::string line;
    auto* coordinates(new double[dim]);
    nodes.reserve(n_nodes);

    std::size_t k(0);
    while (k < n_nodes && !ins.fail())
    {
        getline(ins, line);
        if (ins.fail())
        {
            ERR("TetGenInterface::parseNodes(): Error reading node %d.", k);
            return false;
        }

        std::size_t id;
        std::size_t pos_end = 0;
        std::size_t pos_beg = line.find_first_not_of(' ', pos_end);
        pos_end = line.find_first_of(" \n", pos_beg);

        if (line.empty() || pos_beg==pos_end || line.compare(pos_beg,1,"#") == 0)
        {
            continue;
        }

        if (pos_beg != std::string::npos && pos_end != std::string::npos) {
            id = BaseLib::str2number<std::size_t> (line.substr(pos_beg, pos_end - pos_beg));
            if (k == 0 && id == 0)
                _zero_based_idx = true;
        } else {
            ERR("TetGenInterface::parseNodes(): Error reading ID of node %d.", k);
            delete [] coordinates;
            return false;
        }
        // read coordinates
        const unsigned offset = (_zero_based_idx) ? 0 : 1;
        for (std::size_t i(0); i < dim; i++) {
            pos_beg = line.find_first_not_of(' ', pos_end);
            pos_end = line.find_first_of(" \n", pos_beg);
            if (pos_end == std::string::npos) pos_end = line.size();
            if (pos_beg != std::string::npos)
                coordinates[i] = BaseLib::str2number<double>(line.substr(pos_beg, pos_end-pos_beg));
            else {
                ERR("TetGenInterface::parseNodes(): error reading coordinate %d of node %d.", i, k);
                delete [] coordinates;
                return false;
            }
        }

        nodes.push_back(new MeshLib::Node(coordinates, id-offset));
        // read attributes and boundary markers ... - at the moment we do not use this information
        ++k;
    }

    delete [] coordinates;
    return true;
}

bool TetGenInterface::readElementsFromStream(std::ifstream &ins,
                                             std::vector<MeshLib::Element*> &elements,
                                             std::vector<int> &materials,
                                             const std::vector<MeshLib::Node*> &nodes) const
{
    std::string line;
    getline (ins, line);
    std::size_t n_tets, n_nodes_per_tet;
    bool region_attributes;

    while (!ins.fail())
    {
        BaseLib::simplify(line);
        if (line.empty() || line.compare(0,1,"#") == 0)
        {
            // this line is a comment - skip
            getline (ins, line);
            continue;
        }

        // read header line
        bool header_okay = parseElementsFileHeader(line, n_tets, n_nodes_per_tet, region_attributes);
        if (!header_okay)
            return false;
        if (!parseElements(ins, elements, materials, nodes, n_tets, n_nodes_per_tet, region_attributes))
            return false;
        return true;
    }
    return false;
}

bool TetGenInterface::parseElementsFileHeader(std::string &line,
                                              std::size_t& n_tets,
                                              std::size_t& n_nodes_per_tet,
                                              bool& region_attribute) const
{
    std::size_t pos_beg, pos_end;

    // number of tetrahedras
    pos_beg = line.find_first_not_of (' ');
    pos_end = line.find_first_of(' ', pos_beg);
    if (pos_beg != std::string::npos && pos_end != std::string::npos)
        n_tets = BaseLib::str2number<std::size_t> (line.substr(pos_beg, pos_end - pos_beg));
    else {
        ERR("TetGenInterface::parseElementsFileHeader(): Could not read number of tetrahedra specified in header.");
        return false;
    }
    // nodes per tet - either 4 or 10
    pos_beg = line.find_first_not_of (" \t", pos_end);
    pos_end = line.find_first_of(" \t", pos_beg);
    n_nodes_per_tet = BaseLib::str2number<std::size_t> (line.substr(pos_beg, pos_end - pos_beg));
    // region attribute at tetrahedra?
    pos_beg = line.find_first_not_of (" \t", pos_end);
    pos_end = line.find_first_of(" \t\n", pos_beg);
    if (pos_end == std::string::npos)
        pos_end = line.size();
    region_attribute = line.substr(pos_beg, pos_end - pos_beg) == "1";

    return true;
}

bool TetGenInterface::parseElements(std::ifstream& ins,
                                    std::vector<MeshLib::Element*> &elements,
                                    std::vector<int> &materials,
                                    const std::vector<MeshLib::Node*> &nodes,
                                    std::size_t n_tets,
                                    std::size_t n_nodes_per_tet,
                                    bool region_attribute) const
{
    std::string line;
    auto* ids(static_cast<std::size_t*>(
        alloca(sizeof(std::size_t) * n_nodes_per_tet)));
    elements.reserve(n_tets);
    materials.reserve(n_tets);

    const unsigned offset = (_zero_based_idx) ? 0 : 1;
    for (std::size_t k(0); k < n_tets && !ins.fail(); k++)
    {
        getline (ins, line);
        if (ins.fail())
        {
            ERR("TetGenInterface::parseElements(): Error reading tetrahedron %d.", k);
            return false;
        }

        std::size_t pos_end = 0;
        std::size_t pos_beg = line.find_first_not_of(' ', pos_end);
        pos_end = line.find_first_of(" \n", pos_beg);

        if (line.empty() || pos_beg==pos_end || line.compare(pos_beg,1,"#") == 0)
        {
            k--;
            continue;
        }

        if (pos_beg == std::string::npos || pos_end == std::string::npos)
        {
            ERR("TetGenInterface::parseElements(): Error reading id of tetrahedron %d.", k);
            return false;
        }

        // read node ids
        for (std::size_t i(0); i < n_nodes_per_tet; i++)
        {
            pos_beg = line.find_first_not_of(' ', pos_end);
            pos_end = line.find_first_of(' ', pos_beg);
            if (pos_end == std::string::npos)
                pos_end = line.size();
            if (pos_beg != std::string::npos && pos_end != std::string::npos)
                ids[i] = BaseLib::str2number<std::size_t>(line.substr(pos_beg, pos_end - pos_beg)) - offset;
            else
            {
                ERR("TetGenInterface::parseElements(): Error reading node %d of tetrahedron %d.", i, k);
                return false;
            }
        }

        // read region attribute - this is something like material group
        int region (0);
        if (region_attribute) {
            pos_beg = line.find_first_not_of(' ', pos_end);
            pos_end = line.find_first_of(' ', pos_beg);
            if (pos_end == std::string::npos) pos_end = line.size();
            if (pos_beg != std::string::npos && pos_end != std::string::npos)
                region = BaseLib::str2number<int> (line.substr(pos_beg, pos_end - pos_beg));
            else {
                ERR("TetGenInterface::parseElements(): Error reading region attribute of tetrahedron %d.", k);
                return false;
            }
        }
        // insert new element into vector
        auto** tet_nodes = new MeshLib::Node*[4];
        for (unsigned k(0); k<4; k++) {
            tet_nodes[k] = nodes[ids[k]];
        }
        elements.push_back (new MeshLib::Tet(tet_nodes));
        materials.push_back(region);
    }

    return true;
}

bool TetGenInterface::writeTetGenSmesh(const std::string &file_name,
                                       const GeoLib::GEOObjects &geo_objects,
                                       const std::string &geo_name,
                                       const std::vector<GeoLib::Point> &attribute_points) const
{
    std::vector<GeoLib::Point*> const*const points = geo_objects.getPointVec(geo_name);
    std::vector<GeoLib::Surface*> const*const surfaces = geo_objects.getSurfaceVec(geo_name);

    if (points==nullptr)
    {
        ERR ("Geometry %s not found.", geo_name.c_str());
        return false;
    }
    if (surfaces==nullptr)
        WARN ("No surfaces found for geometry %s. Writing points only.", geo_name.c_str());

    std::ofstream out( file_name.c_str(), std::ios::out );
    out.precision(std::numeric_limits<double>::digits10);
    // the points header
    const std::size_t nPoints (points->size());
    out << nPoints << " 3\n";
    // the point list
    for (std::size_t i=0; i<nPoints; ++i)
        out << i << "  " << (*(*points)[i])[0] << " " << (*(*points)[i])[1] << " " << (*(*points)[i])[2] << "\n";
    // the surfaces header
    const std::size_t nSurfaces = (surfaces) ? surfaces->size() : 0;
    std::size_t nTotalTriangles (0);
    for (std::size_t i=0; i<nSurfaces; ++i)
        nTotalTriangles += (*surfaces)[i]->getNumberOfTriangles();
    out << nTotalTriangles << " 1\n";

    for (std::size_t i=0; i<nSurfaces; ++i)
    {
        const std::size_t nTriangles ((*surfaces)[i]->getNumberOfTriangles());
        const std::size_t marker (i+1); // must NOT be 0!
        // the poly list
        for (std::size_t j=0; j<nTriangles; ++j)
        {
            const GeoLib::Triangle &tri = *(*(*surfaces)[i])[j];
            out << "3  " << tri[0] << " " << tri[1] << " " << tri[2] << " " << marker << "\n";
        }
    }
    out << "0\n"; // the polygon holes list
    // the region attributes list
    if (attribute_points.empty())
        out << "0\n";
    else
    {
        const std::size_t nAttributePoints (attribute_points.size());
        out << nAttributePoints << "\n";
        for (std::size_t i=0; i<nAttributePoints; ++i)
            out << i+1 << " " << attribute_points[i][0] << " " << attribute_points[i][1] << " " << attribute_points[i][2] << " " << 10*attribute_points[i].getID() << "\n";
    }
    INFO ("TetGenInterface::writeTetGenSmesh() - %d points and %d surfaces successfully written.", nPoints, nSurfaces);
    out.close();
    return true;
}

bool TetGenInterface::writeTetGenSmesh(const std::string &file_name,
                                       const MeshLib::Mesh &mesh,
                                       std::vector<MeshLib::Node> &attribute_points) const
{
    if (mesh.getDimension() == 1)
        return false;

    const std::vector<MeshLib::Node*> &nodes = mesh.getNodes();

    std::ofstream out( file_name.c_str(), std::ios::out );
    out.precision(std::numeric_limits<double>::digits10);
    // the points header
    const std::size_t nPoints (nodes.size());
    out << nPoints << " 3\n";
    // the point list
    for (std::size_t i=0; i<nPoints; ++i)
        out << i << "  " << (*nodes[i])[0] << " " << (*nodes[i])[1] << " " << (*nodes[i])[2] << "\n";

    if (mesh.getDimension() == 2)
        write2dElements(out, mesh);
    else
        write3dElements(out, mesh, attribute_points);

    out << "0\n"; // the polygon holes list

    // the region attributes list
    if (attribute_points.empty())
        out << "0\n";
    else
    {
        const std::size_t nAttributePoints (attribute_points.size());
        out << nAttributePoints << "\n";
        for (std::size_t i=0; i<nAttributePoints; ++i)
            out << i+1 << " " << attribute_points[i][0] << " " << attribute_points[i][1] << " " << attribute_points[i][2] << " " << 10*attribute_points[i].getID() << "\n";
    }

    INFO ("TetGenInterface::writeTetGenPoly() - %d points and %d surfaces successfully written.", nPoints, mesh.getNumberOfElements());
    out.close();
    return true;
}

void TetGenInterface::write2dElements(std::ofstream &out,
                                      const MeshLib::Mesh &mesh) const
{
    // the surfaces header
    const std::array<unsigned,7> types = MeshLib::MeshInformation::getNumberOfElementTypes(mesh);
    const unsigned nTotalTriangles (types[1] + (2*types[2]));
    out << nTotalTriangles << " 1\n";

    const std::vector<MeshLib::Element*> &elements = mesh.getElements();
    MeshLib::PropertyVector<int> const*const mat_ids =
        mesh.getProperties().existsPropertyVector<int>("MaterialIDs")
            ? mesh.getProperties().getPropertyVector<int>("MaterialIDs")
            : nullptr;
    const std::size_t nElements (elements.size());
    unsigned element_count(0);
    for (std::size_t i=0; i<nElements; ++i)
    {
        std::string const matId = mat_ids ? std::to_string((*mat_ids)[i]) : "";
        this->writeElementToFacets(out, *elements[i], element_count, matId);
    }
}

void TetGenInterface::write3dElements(std::ofstream &out,
                                      const MeshLib::Mesh &mesh,
                                      std::vector<MeshLib::Node> &attribute_points) const
{
    const std::vector<MeshLib::Element*> &elements = mesh.getElements();
    const std::size_t nElements (elements.size());
    if (!attribute_points.empty())
        attribute_points.clear();

    // get position where number of facets need to be written and figure out worst case of chars that are needed
    const std::streamoff before_elems_pos (out.tellp());
    const unsigned n_spaces (static_cast<unsigned>(std::floor(log(nElements*8))) + 1);
    out << std::string(n_spaces, ' ') << " 1\n";
    auto const* const materialIds =
        mesh.getProperties().getPropertyVector<int>("MaterialIDs");
    unsigned element_count(0);
    for (std::size_t i=0; i<nElements; ++i)
    {
        if (elements[i]->getDimension() < 3)
            continue;

        const unsigned nFaces (elements[i]->getNumberOfNeighbors());
        std::string const mat_id_str =
            materialIds ? std::to_string((*materialIds)[i]) : "";
        for (std::size_t j=0; j<nFaces; ++j)
        {
            MeshLib::Element const*const neighbor ( elements[i]->getNeighbor(j) );

            if (neighbor && materialIds && (*materialIds)[i] <= (*materialIds)[neighbor->getID()])
                continue;

            std::unique_ptr<MeshLib::Element const> const face (elements[i]->getFace(j));
            this->writeElementToFacets(out, *face, element_count, mat_id_str);
        }
        if (materialIds)
            attribute_points.emplace_back(
                elements[i]->getCenterOfGravity().getCoords(),
                (*materialIds)[i]);
    }
    // add number of facets at correct position and jump back
    const std::streamoff after_elems_pos (out.tellp());
    out.seekp(before_elems_pos);
    out << element_count;
    out.seekp(after_elems_pos);
}

void TetGenInterface::writeElementToFacets(std::ofstream &out, const MeshLib::Element &element, unsigned &element_count, std::string const& matId) const
{
    element_count++;
    if (element.getGeomType() == MeshLib::MeshElemType::TRIANGLE)
        out << "3  " << element.getNodeIndex(0) << " " << element.getNodeIndex(1) << " " << element.getNodeIndex(2) << " " << matId << " # " << element_count << "\n";
    else if (element.getGeomType() == MeshLib::MeshElemType::QUAD)
    {
        out << "3  " << element.getNodeIndex(0) << " " << element.getNodeIndex(1) << " " << element.getNodeIndex(2) << " " << matId << " # " << element_count << "\n";
        element_count++;
        out << "3  " << element.getNodeIndex(0) << " " << element.getNodeIndex(2) << " " << element.getNodeIndex(3) << " " << matId << " # " << element_count << "\n";
    }
}

} // end namespace FileIO
