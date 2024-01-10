/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-08
 * \brief  Implementation of the MeshIO class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshIO.h"

#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>

#include "BaseLib/Error.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/StringTools.h"
#include "MeshLib/Elements/Elements.h"
#include "MeshLib/Location.h"
#include "MeshLib/Node.h"
#include "MeshLib/PropertyVector.h"

namespace
{
std::size_t readMaterialID(std::istream& in)
{
    unsigned index;
    unsigned material_id;
    if (!(in >> index >> material_id))
    {
        return std::numeric_limits<std::size_t>::max();
    }
    return material_id;
}

MeshLib::Element* readElement(std::istream& in,
                              const std::vector<MeshLib::Node*>& nodes)
{
    std::string elem_type_str;
    MeshLib::MeshElemType elem_type(MeshLib::MeshElemType::INVALID);

    do
    {
        if (!(in >> elem_type_str))
        {
            return nullptr;
        }
        elem_type = MeshLib::String2MeshElemType(elem_type_str);
    } while (elem_type == MeshLib::MeshElemType::INVALID);

    auto* idx = new unsigned[8];
    MeshLib::Element* elem;

    switch (elem_type)
    {
        case MeshLib::MeshElemType::LINE:
        {
            for (int i = 0; i < 2; ++i)
            {
                if (!(in >> idx[i]))
                {
                    return nullptr;
                }
            }
            // edge_nodes array will be deleted from Line object
            auto** edge_nodes = new MeshLib::Node*[2];
            for (unsigned k(0); k < 2; ++k)
            {
                edge_nodes[k] = nodes[idx[k]];
            }
            elem = new MeshLib::Line(edge_nodes);
            break;
        }
        case MeshLib::MeshElemType::TRIANGLE:
        {
            for (int i = 0; i < 3; ++i)
            {
                if (!(in >> idx[i]))
                {
                    return nullptr;
                }
            }
            auto** tri_nodes = new MeshLib::Node*[3];
            for (unsigned k(0); k < 3; ++k)
            {
                tri_nodes[k] = nodes[idx[k]];
            }
            elem = new MeshLib::Tri(tri_nodes);
            break;
        }
        case MeshLib::MeshElemType::QUAD:
        {
            for (int i = 0; i < 4; ++i)
            {
                if (!(in >> idx[i]))
                {
                    return nullptr;
                }
            }
            auto** quad_nodes = new MeshLib::Node*[4];
            for (unsigned k(0); k < 4; ++k)
            {
                quad_nodes[k] = nodes[idx[k]];
            }
            elem = new MeshLib::Quad(quad_nodes);
            break;
        }
        case MeshLib::MeshElemType::TETRAHEDRON:
        {
            for (int i = 0; i < 4; ++i)
            {
                if (!(in >> idx[i]))
                {
                    return nullptr;
                }
            }
            auto** tet_nodes = new MeshLib::Node*[4];
            for (unsigned k(0); k < 4; ++k)
            {
                tet_nodes[k] = nodes[idx[k]];
            }
            elem = new MeshLib::Tet(tet_nodes);
            break;
        }
        case MeshLib::MeshElemType::HEXAHEDRON:
        {
            for (int i = 0; i < 8; ++i)
            {
                if (!(in >> idx[i]))
                {
                    return nullptr;
                }
            }
            auto** hex_nodes = new MeshLib::Node*[8];
            for (unsigned k(0); k < 8; ++k)
            {
                hex_nodes[k] = nodes[idx[k]];
            }
            elem = new MeshLib::Hex(hex_nodes);
            break;
        }
        case MeshLib::MeshElemType::PYRAMID:
        {
            for (int i = 0; i < 5; ++i)
            {
                if (!(in >> idx[i]))
                {
                    return nullptr;
                }
            }
            auto** pyramid_nodes = new MeshLib::Node*[5];
            for (unsigned k(0); k < 5; ++k)
            {
                pyramid_nodes[k] = nodes[idx[k]];
            }
            elem = new MeshLib::Pyramid(pyramid_nodes);
            break;
        }
        case MeshLib::MeshElemType::PRISM:
        {
            for (int i = 0; i < 6; ++i)
            {
                if (!(in >> idx[i]))
                {
                    return nullptr;
                }
            }
            auto** prism_nodes = new MeshLib::Node*[6];
            for (unsigned k(0); k < 6; ++k)
            {
                prism_nodes[k] = nodes[idx[k]];
            }
            elem = new MeshLib::Prism(prism_nodes);
            break;
        }
        default:
            elem = nullptr;
            break;
    }

    delete[] idx;

    return elem;
}

std::string ElemType2StringOutput(const MeshLib::MeshElemType t)
{
    if (t == MeshLib::MeshElemType::LINE)
    {
        return "line";
    }
    if (t == MeshLib::MeshElemType::QUAD)
    {
        return "quad";
    }
    if (t == MeshLib::MeshElemType::HEXAHEDRON)
    {
        return "hex";
    }
    if (t == MeshLib::MeshElemType::TRIANGLE)
    {
        return "tri";
    }
    if (t == MeshLib::MeshElemType::TETRAHEDRON)
    {
        return "tet";
    }
    if (t == MeshLib::MeshElemType::PRISM)
    {
        return "pris";
    }
    if (t == MeshLib::MeshElemType::PYRAMID)
    {
        return "pyra";
    }
    return "none";
}

void writeElements(std::vector<MeshLib::Element*> const& ele_vec,
                   MeshLib::PropertyVector<int> const* const material_ids,
                   std::ostream& out)
{
    const std::size_t ele_vector_size(ele_vec.size());

    out << ele_vector_size << "\n";
    for (std::size_t i(0); i < ele_vector_size; ++i)
    {
        auto const& element = *ele_vec[i];
        if (element.getNumberOfBaseNodes() != element.getNumberOfNodes())
        {
            OGS_FATAL(
                "Found high order element in the mesh that is not required by "
                "OGS 5 input. High order elements are generated in OGS 5 on "
                "demand.");
        }

        out << i << " ";
        if (!material_ids)
        {
            out << "0 ";
        }
        else
        {
            out << (*material_ids)[i] << " ";
        }
        out << ElemType2StringOutput(element.getGeomType()) << " ";
        unsigned nElemNodes(element.getNumberOfBaseNodes());
        for (std::size_t j = 0; j < nElemNodes; ++j)
        {
            out << element.getNode(j)->getID() << " ";
        }
        out << "\n";
    }
}

}  // namespace

namespace MeshLib
{
namespace IO
{
namespace Legacy
{
MeshIO::MeshIO() = default;

MeshLib::Mesh* MeshIO::loadMeshFromFile(const std::string& file_name)
{
    INFO("Reading OGS legacy mesh ... ");

    std::ifstream in(file_name.c_str(), std::ios::in);
    if (!in.is_open())
    {
        WARN("MeshIO::loadMeshFromFile() - Could not open file {:s}.",
             file_name);
        return nullptr;
    }

    std::string line_string;
    std::getline(in, line_string);

    if (line_string.find("#FEM_MSH") != std::string::npos)  // OGS mesh file
    {
        std::vector<MeshLib::Node*> nodes;
        std::vector<MeshLib::Element*> elements;
        std::vector<std::size_t> materials;

        while (!in.eof())
        {
            std::getline(in, line_string);

            // check keywords
            if (line_string.find("#STOP") != std::string::npos)
            {
                break;
            }
            if (line_string.find("$NODES") != std::string::npos)
            {
                std::getline(in, line_string);
                BaseLib::trim(line_string);
                unsigned nNodes = atoi(line_string.c_str());
                for (unsigned i = 0; i < nNodes; ++i)
                {
                    std::getline(in, line_string);
                    std::stringstream iss(line_string);
                    unsigned idx;
                    double x;
                    double y;
                    double z;
                    iss >> idx >> x >> y >> z;
                    auto* node(new MeshLib::Node(x, y, z, idx));
                    nodes.push_back(node);
                    std::string s;
                    iss >> s;
                    if (s.find("$AREA") != std::string::npos)
                    {
                        double double_dummy;
                        iss >> double_dummy;
                    }
                }
            }
            else if (line_string.find("$ELEMENTS") != std::string::npos)
            {
                std::getline(in, line_string);
                BaseLib::trim(line_string);
                unsigned nElements = atoi(line_string.c_str());
                for (unsigned i = 0; i < nElements; ++i)
                {
                    std::getline(in, line_string);
                    std::stringstream ss(line_string);
                    materials.push_back(readMaterialID(ss));
                    MeshLib::Element* elem(readElement(ss, nodes));
                    if (elem == nullptr)
                    {
                        ERR("Reading mesh element {:d} from file '{:s}' "
                            "failed.",
                            i, file_name);
                        // clean up the elements vector
                        std::for_each(elements.begin(), elements.end(),
                                      std::default_delete<MeshLib::Element>());
                        // clean up the nodes vector
                        std::for_each(nodes.begin(), nodes.end(),
                                      std::default_delete<MeshLib::Node>());
                        return nullptr;
                    }
                    elements.push_back(elem);
                }
            }
        }

        if (elements.empty())
        {
            ERR("MeshIO::loadMeshFromFile() - File did not contain element "
                "information.");
            for (auto& node : nodes)
            {
                delete node;
            }
            return nullptr;
        }

        MeshLib::Mesh* mesh(new MeshLib::Mesh(
            BaseLib::extractBaseNameWithoutExtension(file_name), nodes,
            elements, true /* compute_element_neighbors */));

        auto* const material_ids =
            mesh->getProperties().createNewPropertyVector<int>(
                "MaterialIDs", MeshLib::MeshItemType::Cell, 1);
        if (!material_ids)
        {
            WARN("Could not create PropertyVector for MaterialIDs in Mesh.");
        }
        else
        {
            material_ids->insert(material_ids->end(), materials.cbegin(),
                                 materials.cend());
        }
        INFO("\t... finished.");
        INFO("Nr. Nodes: {:d}.", nodes.size());
        INFO("Nr. Elements: {:d}.", elements.size());

        in.close();
        return mesh;
    }

    in.close();
    return nullptr;
}

bool MeshIO::write()
{
    if (!_mesh)
    {
        WARN("MeshIO::write(): Cannot write: no mesh object specified.");
        return false;
    }

    out << "#FEM_MSH\n"
        << "$PCS_TYPE\n"
        << "  NO_PCS\n"
        << "$NODES\n"
        << "  ";
    const std::size_t n_nodes(_mesh->getNumberOfNodes());
    out << n_nodes << "\n";
    for (std::size_t i(0); i < n_nodes; ++i)
    {
        out << i << " " << *(_mesh->getNode(i)) << "\n";
    }

    out << "$ELEMENTS\n"
        << "  ";

    if (!_mesh->getProperties().existsPropertyVector<int>("MaterialIDs"))
    {
        writeElements(_mesh->getElements(), nullptr, out);
    }
    else
    {
        writeElements(
            _mesh->getElements(),
            _mesh->getProperties().getPropertyVector<int>("MaterialIDs"), out);
    }
    out << "#STOP\n";

    return true;
}

void MeshIO::setMesh(const MeshLib::Mesh* mesh)
{
    _mesh = mesh;
}

}  // end namespace Legacy
}  // end namespace IO
}  // namespace MeshLib
