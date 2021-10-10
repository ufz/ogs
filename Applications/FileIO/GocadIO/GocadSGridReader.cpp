/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "GocadSGridReader.h"

#include <algorithm>
#include <boost/tokenizer.hpp>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>

#include "BaseLib/FileTools.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Mesh.h"

namespace FileIO
{
namespace Gocad
{
using Bitset = boost::dynamic_bitset<>;

GocadSGridReader::GocadSGridReader(std::string const& fname)
    : _fname(fname),
      _path(BaseLib::extractPath(fname)),
      _n_face_sets(0),
      _double_precision_binary(false),
      _bin_pnts_in_double_precision(false)
{
    // check if file exists
    std::ifstream in(_fname.c_str());
    if (!in)
    {
        ERR("Could not open '{:s}'.", _fname);
        in.close();
        return;
    }

    bool pnts_read(false);

    // read information in the stratigraphic grid file
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty())  // skip empty lines
        {
            continue;
        }
        if (line.back() == '\r')  // check dos line ending
        {
            line.pop_back();
        }
        if (line.substr(0, 8) == "HEADER {")
        {
            parseHeader(in);
        }
        if (line.substr(0, 32) == "GOCAD_ORIGINAL_COORDINATE_SYSTEM")
        {
            _coordinate_system.parse(in);
            continue;
        }
        if (line.substr(0, 7) == "AXIS_N ")
        {
            parseDims(line);
        }
        else if (line.substr(0, 12) == "POINTS_FILE ")
        {
            parseFileName(line, _pnts_fname);
        }
        else if (line.substr(0, 9) == "PROPERTY ")
        {
            _property_meta_data_vecs.push_back(
                Gocad::parseGocadPropertyMetaData(line, in, _path));
        }
        else if (line.substr(0, 35) == "BINARY_POINTS_IN_DOUBLE_PRECISION 1")
        {
            _bin_pnts_in_double_precision = true;
        }
        else if (line.substr(0, 11) == "FLAGS_FILE ")
        {
            parseFileName(line, _flags_fname);
        }
        else if (line.substr(0, 18) == "REGION_FLAGS_FILE ")
        {
            parseFileName(line, _region_flags_fname);
        }
        else if (line.substr(0, 7) == "REGION " ||
                 line.substr(0, 13) == "MODEL_REGION ")
        {
            regions.push_back(Gocad::parseRegion(line));
        }
        else if (line.substr(0, 12) == "MODEL_LAYER ")
        {
            layers.push_back(Gocad::parseLayer(line, regions));
        }
        else if (line.substr(0, 24) == "REGION_FLAGS_BIT_LENGTH ")
        {
            std::istringstream iss(line);
            std::istream_iterator<std::string> it(iss);
            it++;
            std::size_t bit_length = std::atoi(it->c_str());
            if (regions.size() != bit_length)
            {
                ERR("{:d} regions read but {:d} expected.\n", regions.size(),
                    bit_length);
                throw std::runtime_error(
                    "Number of read regions differs from expected.\n");
            }
        }
        else if (line.substr(0, 9) == "FACE_SET ")
        {
            // first read the points
            if (!pnts_read)
            {
                readNodesBinary();
                pnts_read = true;
            }
            parseFaceSet(line, in);
        }
    }

#ifndef NDEBUG
    std::stringstream regions_ss;
    regions_ss << regions.size() << " regions read:\n";
    std::copy(regions.begin(), regions.end(),
              std::ostream_iterator<Gocad::Region>(regions_ss, "\t"));
    DBUG("{:s}", regions_ss.str());

    std::stringstream layers_ss;
    layers_ss << layers.size() << " layers read:\n";
    std::copy(layers.begin(), layers.end(),
              std::ostream_iterator<Gocad::Layer>(layers_ss, "\n"));
    DBUG("{:s}", layers_ss.str());

    std::stringstream properties_ss;
    properties_ss << "meta data for " << _property_meta_data_vecs.size()
                  << " properties read:\n";
    std::copy(_property_meta_data_vecs.begin(), _property_meta_data_vecs.end(),
              std::ostream_iterator<Gocad::Property>(properties_ss, "\n"));
    DBUG("{:s}", properties_ss.str());
#endif

    // if not done already read the points
    if (!pnts_read)
    {
        readNodesBinary();
    }
    readElementPropertiesBinary();
    std::vector<Bitset> region_flags = readRegionFlagsBinary();
    mapRegionFlagsToCellProperties(region_flags);

    readSplitInformation();

    in.close();
}

GocadSGridReader::~GocadSGridReader()
{
    for (auto node : _nodes)
    {
        delete node;
    }
    for (auto node : _split_nodes)
    {
        delete node;
    }
}

std::unique_ptr<MeshLib::Mesh> GocadSGridReader::getMesh() const
{
    std::vector<MeshLib::Node*> nodes;
    std::transform(_nodes.cbegin(), _nodes.cend(), std::back_inserter(nodes),
                   [](MeshLib::Node const* const node)
                   { return new MeshLib::Node(*node); });

    std::vector<MeshLib::Element*> elements(createElements(nodes));
    applySplitInformation(nodes, elements);

    DBUG("Creating mesh from Gocad SGrid.");
    std::unique_ptr<MeshLib::Mesh> mesh(new MeshLib::Mesh(
        BaseLib::extractBaseNameWithoutExtension(_fname), nodes, elements));
    addGocadPropertiesToMesh(*mesh);
    DBUG("Mesh created.");

    return mesh;
}

void GocadSGridReader::addGocadPropertiesToMesh(MeshLib::Mesh& mesh) const
{
    std::vector<std::string> const& prop_names(getPropertyNames());
    for (auto const& name : prop_names)
    {
        auto prop = getProperty(name);
        if (!prop)
        {
            continue;
        }

        DBUG("Adding Gocad property '{:s}' with {:d} values.", name,
             prop->_property_data.size());

        auto pv = MeshLib::getOrCreateMeshProperty<double>(
            mesh, name, MeshLib::MeshItemType::Cell, 1);
        if (pv == nullptr)
        {
            ERR("Could not create mesh property '{:s}'.", name);
            continue;
        }

        pv->resize(prop->_property_data.size());
        std::copy(prop->_property_data.cbegin(), prop->_property_data.cend(),
                  pv->begin());
    }
}

void GocadSGridReader::parseHeader(std::istream& in)
{
    std::string line;
    while (std::getline(in, line))
    {
        if (line.front() == '}')
        {
            return;
        }
        if (line.substr(0, 27) == "double_precision_binary: on")
        {
            _double_precision_binary = true;
        }
    }
    if (_double_precision_binary)
    {
        DBUG(
            "GocadSGridReader::parseHeader(): _double_precision_binary == "
            "true.");
    }
}

void GocadSGridReader::parseDims(std::string const& line)
{
    std::size_t x_dim(0);
    std::size_t y_dim(0);
    std::size_t z_dim(0);
    boost::tokenizer<> tok(line);
    auto it(tok.begin());
    it++;  // overread token "AXIS"
    it++;  // overread "N"
    std::stringstream ssx(*(it),
                          std::stringstream::in | std::stringstream::out);
    ssx >> x_dim;
    it++;
    std::stringstream ssy(*it, std::stringstream::in | std::stringstream::out);
    ssy >> y_dim;
    it++;
    std::stringstream ssz(*it, std::stringstream::in | std::stringstream::out);
    ssz >> z_dim;
    _index_calculator = Gocad::IndexCalculator(x_dim, y_dim, z_dim);
    DBUG(
        "x_dim = {:d}, y_dim = {:d}, z_dim = {:d} => #nodes = {:d}, #cells = "
        "{:d}",
        x_dim, y_dim, z_dim, _index_calculator._n_nodes,
        _index_calculator._n_cells);
}

void GocadSGridReader::parseFileName(std::string const& line,
                                     std::string& result_string) const
{
    boost::char_separator<char> sep(" \r");
    boost::tokenizer<boost::char_separator<char>> tok(line, sep);
    auto it(tok.begin());
    ++it;  // overread POINTS_FILE or FLAGS_FILE or REGION_FLAGS_FILE
    result_string = _path + *it;
}

/**
 * @param line input/output
 * @param in input stream containing the face set
 */
void GocadSGridReader::parseFaceSet(std::string& line, std::istream& in)
{
    // create and initialize a Gocad::Property object for storing face set data
    Gocad::Property face_set_property;
    face_set_property._property_id = _n_face_sets;
    face_set_property._property_name = "FaceSet";
    face_set_property._property_class_name = "FaceSetData";
    face_set_property._property_unit = "unitless";
    face_set_property._property_data_type = "double";
    face_set_property._property_data_fname = "";
    face_set_property._property_no_data_value = -1.0;
    face_set_property._property_data.resize(_index_calculator._n_cells);
    std::fill(face_set_property._property_data.begin(),
              face_set_property._property_data.end(),
              face_set_property._property_no_data_value);

    std::istringstream iss(line);
    std::istream_iterator<std::string> it(iss);
    // Check first word is FACE_SET
    if (*it != std::string("FACE_SET"))
    {
        ERR("Expected FACE_SET keyword but '{:s}' found.", it->c_str());
        throw std::runtime_error(
            "In GocadSGridReader::parseFaceSet() expected FACE_SET keyword not "
            "found.");
    }
    ++it;
    face_set_property._property_name += *it;
    ++it;
    auto const n_of_face_set_ids(
        static_cast<std::size_t>(std::atoi(it->c_str())));
    std::size_t face_set_id_cnt(0);

    while (std::getline(in, line) && face_set_id_cnt < n_of_face_set_ids)
    {
        if (line.back() == '\r')
        {
            line.pop_back();
        }
        boost::char_separator<char> sep("\t ");
        boost::tokenizer<boost::char_separator<char>> tokens(line, sep);

        for (auto tok_it = tokens.begin(); tok_it != tokens.end();)
        {
            auto const id(static_cast<std::size_t>(std::atoi(tok_it->c_str())));
            tok_it++;
            auto const face_indicator(
                static_cast<std::size_t>(std::atoi(tok_it->c_str())));
            tok_it++;

            if (id >= _index_calculator._n_nodes)
            {
                ERR("Face set id {:d} is greater than the number of nodes "
                    "({:d}).",
                    id, _index_calculator._n_nodes);
            }
            else
            {
                static_cast<GocadNode*>(_nodes[id])
                    ->setFaceSet(_n_face_sets, face_indicator);
                std::array<std::size_t, 3> const c(
                    _index_calculator.getCoordsForID(id));
                if (c[0] >= _index_calculator._x_dim - 1)
                {
                    ERR("****** i coord {:d} to big for id {:d}.", c[0], id);
                }
                if (c[1] >= _index_calculator._y_dim - 1)
                {
                    ERR("****** j coord {:d} to big for id {:d}.", c[1], id);
                }
                if (c[2] >= _index_calculator._z_dim - 1)
                {
                    ERR("****** k coord {:d} to big for id {:d}.", c[2], id);
                }
                std::size_t const cell_id(
                    _index_calculator.getCellIdx(c[0], c[1], c[2]));
                face_set_property._property_data[cell_id] =
                    static_cast<double>(face_indicator);
            }
            face_set_id_cnt++;
        }
    }

    if (face_set_id_cnt != n_of_face_set_ids)
    {
        ERR("Expected {:d} number of face set ids, read {:d}.",
            n_of_face_set_ids, face_set_id_cnt);
        throw std::runtime_error(
            "Expected number of face set points does not match number of read "
            "points.");
    }
    _n_face_sets++;

    // pre condition: split nodes are read already
    for (auto split_node : _split_nodes)
    {
        std::size_t const id(_index_calculator(split_node->getGridCoords()));
        split_node->transmitFaceDirections(*_nodes[id]);
    }

    _property_meta_data_vecs.push_back(face_set_property);
}

// Reads given number of bits (rounded up to next byte) into a bitset.
// Used for reading region information which can be represented by some
// number of bits.
Bitset readBits(std::ifstream& in, const std::size_t bits)
{
    using block_t = Bitset::block_type;
    auto const bytes = static_cast<std::size_t>(std::ceil(bits / 8.));
    std::size_t const blocks =
        bytes + 1 < sizeof(block_t) ? 1 : (bytes + 1) / sizeof(block_t);

    std::vector<block_t> data;
    data.resize(blocks);
    std::fill_n(data.data(), blocks, 0);
    in.read(reinterpret_cast<char*>(data.data()), bytes);

    return Bitset(data.begin(), data.end());
}

void GocadSGridReader::readNodesBinary()
{
    std::ifstream in(_pnts_fname.c_str(), std::ios::binary);
    if (!in)
    {
        ERR("Could not open points file '{:s}'.", _pnts_fname);
        throw std::runtime_error("Could not open points file.");
    }

    std::size_t const n = _index_calculator._n_nodes;
    _nodes.resize(n);

    double coords[3];

    std::size_t k = 0;
    while (in && k < n * 3)
    {
        if (_bin_pnts_in_double_precision)
        {
            coords[k % 3] =
                BaseLib::swapEndianness(BaseLib::readBinaryValue<double>(in));
        }
        else
        {
            coords[k % 3] =
                BaseLib::swapEndianness(BaseLib::readBinaryValue<float>(in));
        }
        if ((k + 1) % 3 == 0)
        {
            const std::size_t layer_transition_idx(
                _index_calculator.getCoordsForID(k / 3)[2]);
            _nodes[k / 3] = new GocadNode(coords, k / 3, layer_transition_idx);
        }
        k++;
    }
    if (k != n * 3 && !in.eof())
    {
        ERR("Read different number of points. Expected {:d} floats, got "
            "{:d}.\n",
            n * 3, k);
    }
}

void GocadSGridReader::mapRegionFlagsToCellProperties(
    std::vector<Bitset> const& rf)
{
    DBUG(
        "GocadSGridReader::mapRegionFlagsToCellProperties region_flags.size: "
        "{:d}",
        rf.size());

    Gocad::Property region_flags;
    region_flags._property_id = 0;
    region_flags._property_name = "RegionFlags";
    region_flags._property_class_name = "RegionFlags";
    region_flags._property_unit = "unitless";
    region_flags._property_data_type = "int";
    region_flags._property_data_fname = "";
    region_flags._property_no_data_value = -1;
    std::size_t const n = _index_calculator._n_cells;
    region_flags._property_data.resize(n);
    std::fill(region_flags._property_data.begin(),
              region_flags._property_data.end(), -1);

    // region flags are stored in each node ijk and give the region index for
    // the ijk-th cell.
    for (std::size_t i(0); i < _index_calculator._x_dim - 1; i++)
    {
        for (std::size_t j(0); j < _index_calculator._y_dim - 1; j++)
        {
            for (std::size_t k(0); k < _index_calculator._z_dim - 1; k++)
            {
                std::size_t const cell_id(
                    _index_calculator.getCellIdx(i, j, k));
                std::size_t const node_id(_index_calculator({i, j, k}));
                for (auto& region : regions)
                {
                    if (rf[node_id].test(region.bit))
                    {
                        region_flags._property_data[cell_id] += region.bit + 1;
                    }
                }
            }
        }
    }

    _property_meta_data_vecs.push_back(region_flags);
}

void GocadSGridReader::readElementPropertiesBinary()
{
    for (auto& property : _property_meta_data_vecs)
    {
        std::string const& fname(property._property_data_fname);
        if (property._property_data_fname.empty())
        {
            WARN("Empty filename for property {:s}.", property._property_name);
            continue;
        }
        std::vector<float> float_properties =
            BaseLib::readBinaryArray<float>(fname, _index_calculator._n_cells);
        DBUG(
            "GocadSGridReader::readElementPropertiesBinary(): Read {:d} float "
            "properties from binary file.",
            _index_calculator._n_cells);

        std::transform(float_properties.cbegin(), float_properties.cend(),
                       float_properties.begin(),
                       [](float const& val)
                       { return BaseLib::swapEndianness(val); });

        property._property_data.resize(float_properties.size());
        std::copy(float_properties.begin(), float_properties.end(),
                  property._property_data.begin());
        if (property._property_data.empty())
        {
            ERR("Reading of element properties file '{:s}' failed.", fname);
        }
    }
}

std::vector<Bitset> GocadSGridReader::readRegionFlagsBinary() const
{
    std::vector<Bitset> result;

    std::ifstream in(_region_flags_fname.c_str());
    if (!in)
    {
        ERR("readRegionFlagsBinary(): Could not open file '{:s}' for input.\n",
            _region_flags_fname);
        in.close();
        return result;
    }

    std::size_t const n = _index_calculator._n_nodes;
    result.resize(n);

    std::size_t k = 0;
    while (in && k < n)
    {
        result[k++] = readBits(in, regions.size());
    }
    if (k != n && !in.eof())
    {
        ERR("Read different number of values. Expected {:d}, got {:d}.\n", n,
            k);
    }
    return result;
}
std::vector<MeshLib::Element*> GocadSGridReader::createElements(
    std::vector<MeshLib::Node*> const& nodes) const
{
    std::vector<MeshLib::Element*> elements;
    elements.resize(_index_calculator._n_cells);
    std::array<MeshLib::Node*, 8> element_nodes{};
    std::size_t cnt(0);
    for (std::size_t k(0); k < _index_calculator._z_dim - 1; k++)
    {
        for (std::size_t j(0); j < _index_calculator._y_dim - 1; j++)
        {
            for (std::size_t i(0); i < _index_calculator._x_dim - 1; i++)
            {
                element_nodes[0] = nodes[_index_calculator({i, j, k})];
                element_nodes[1] = nodes[_index_calculator({i + 1, j, k})];
                element_nodes[2] = nodes[_index_calculator({i + 1, j + 1, k})];
                element_nodes[3] = nodes[_index_calculator({i, j + 1, k})];
                element_nodes[4] = nodes[_index_calculator({i, j, k + 1})];
                element_nodes[5] = nodes[_index_calculator({i + 1, j, k + 1})];
                element_nodes[6] =
                    nodes[_index_calculator({i + 1, j + 1, k + 1})];
                element_nodes[7] = nodes[_index_calculator({i, j + 1, k + 1})];
                elements[cnt] = new MeshLib::Hex(
                    element_nodes, _index_calculator.getCellIdx(i, j, k));
                cnt++;
            }
        }
    }
    return elements;
}

void GocadSGridReader::readSplitInformation()
{
    std::ifstream in(_fname.c_str());
    if (!in)
    {
        ERR("Could not open '{:s}'.", _fname);
        in.close();
        return;
    }

    // read split information from the stratigraphic grid file
    std::string line;
    std::stringstream ss;
    while (std::getline(in, line))
    {
        std::size_t pos(line.find("SPLIT "));
        if (pos != std::string::npos)
        {
            ss << line.substr(pos + 6, line.size() - (pos + 6));
            // read position in grid
            std::array<std::size_t, 3> grid_coords{};
            ss >> grid_coords[0];
            ss >> grid_coords[1];
            ss >> grid_coords[2];
            // read coordinates for the split node
            double coords[3];
            ss >> coords[0];
            ss >> coords[1];
            ss >> coords[2];
            // read the id
            std::size_t id;
            ss >> id;
            // read the affected cells
            std::bitset<8> affected_cells{};
            for (std::size_t ac = 0; ac < 8; ++ac)
            {
                char bit;
                ss >> bit;
                affected_cells[ac] = bit != '0';
            }
            const std::size_t layer_transition_index(
                _nodes[id]->getLayerTransitionIndex());
            _split_nodes.push_back(new GocadSplitNode(coords, id, grid_coords,
                                                      affected_cells,
                                                      layer_transition_index));
        }
    }
}

void GocadSGridReader::applySplitInformation(
    std::vector<MeshLib::Node*>& nodes,
    std::vector<MeshLib::Element*> const& elements) const
{
    for (auto split_node : _split_nodes)
    {
        std::size_t const new_node_pos(nodes.size());
        nodes.push_back(
            new MeshLib::Node(split_node->getCoords(), new_node_pos));

        // get grid coordinates
        std::array<std::size_t, 3> const& gc(split_node->getGridCoords());
        // get affected cells
        auto const& affected_cells(split_node->getAffectedCells());
        // get mesh node to substitute in elements
        MeshLib::Node const* const node2sub(nodes[_index_calculator(gc)]);

        if (affected_cells[0] && gc[0] < _index_calculator._x_dim - 1 &&
            gc[1] < _index_calculator._y_dim - 1 &&
            gc[2] < _index_calculator._z_dim - 1)
        {
            const std::size_t idx(
                _index_calculator.getCellIdx(gc[0], gc[1], gc[2]));
            modifyElement(elements[idx], node2sub, nodes[new_node_pos]);
        }
        if (affected_cells[1] && gc[0] > 0 &&
            gc[1] < _index_calculator._y_dim - 1 &&
            gc[2] < _index_calculator._z_dim - 1)
        {
            const std::size_t idx(
                _index_calculator.getCellIdx(gc[0] - 1, gc[1], gc[2]));
            modifyElement(elements[idx], node2sub, nodes[new_node_pos]);
        }
        if (affected_cells[2] && gc[1] > 0 &&
            gc[0] < _index_calculator._x_dim - 1 &&
            gc[2] < _index_calculator._z_dim - 1)
        {
            const std::size_t idx(
                _index_calculator.getCellIdx(gc[0], gc[1] - 1, gc[2]));
            modifyElement(elements[idx], node2sub, nodes[new_node_pos]);
        }
        if (affected_cells[3] && gc[0] > 0 && gc[1] > 0 &&
            gc[2] < _index_calculator._z_dim - 1)
        {
            const std::size_t idx(
                _index_calculator.getCellIdx(gc[0] - 1, gc[1] - 1, gc[2]));
            modifyElement(elements[idx], node2sub, nodes[new_node_pos]);
        }
        if (affected_cells[4] && gc[2] > 0 &&
            gc[0] < _index_calculator._x_dim - 1 &&
            gc[1] < _index_calculator._y_dim - 1)
        {
            const std::size_t idx(
                _index_calculator.getCellIdx(gc[0], gc[1], gc[2] - 1));
            modifyElement(elements[idx], node2sub, nodes[new_node_pos]);
        }
        if (affected_cells[5] && gc[0] > 0 && gc[2] > 0 &&
            gc[1] < _index_calculator._y_dim - 1)
        {
            const std::size_t idx(
                _index_calculator.getCellIdx(gc[0] - 1, gc[1], gc[2] - 1));
            modifyElement(elements[idx], node2sub, nodes[new_node_pos]);
        }
        if (affected_cells[6] && gc[1] > 0 && gc[2] > 0 &&
            gc[0] < _index_calculator._x_dim - 1)
        {
            const std::size_t idx(
                _index_calculator.getCellIdx(gc[0], gc[1] - 1, gc[2] - 1));
            modifyElement(elements[idx], node2sub, nodes[new_node_pos]);
        }
        if (affected_cells[7] && gc[0] > 0 && gc[1] > 0 && gc[2] > 0)
        {
            const std::size_t idx(
                _index_calculator.getCellIdx(gc[0] - 1, gc[1] - 1, gc[2] - 1));
            modifyElement(elements[idx], node2sub, nodes[new_node_pos]);
        }
    }
}

void GocadSGridReader::modifyElement(MeshLib::Element* hex,
                                     MeshLib::Node const* node2sub,
                                     MeshLib::Node* substitute_node)
{
    // get the node pointers of the cell
    MeshLib::Node* const* hex_nodes(hex->getNodes());
    // search for the position the split node will be set to
    MeshLib::Node* const* node_pos(
        std::find(hex_nodes, hex_nodes + 8, node2sub));
    // set the split node instead of the node2sub
    if (node_pos != hex_nodes + 8)
    {
        const_cast<MeshLib::Node**>(
            hex_nodes)[std::distance(hex_nodes, node_pos)] = substitute_node;
    }
}

Gocad::Property const* GocadSGridReader::getProperty(
    std::string const& name) const
{
    auto const it(std::find_if(_property_meta_data_vecs.begin(),
                               _property_meta_data_vecs.end(),
                               [&name](Gocad::Property const& p)
                               { return p._property_name == name; }));
    if (it == _property_meta_data_vecs.end())
    {
        return nullptr;
    }
    return &*it;
}

std::vector<std::string> GocadSGridReader::getPropertyNames() const
{
    std::vector<std::string> names;
    std::transform(_property_meta_data_vecs.begin(),
                   _property_meta_data_vecs.end(),
                   std::back_inserter(names),
                   [](Gocad::Property const& p) { return p._property_name; });
    return names;
}

std::unique_ptr<MeshLib::Mesh> GocadSGridReader::getFaceSetMesh(
    std::size_t const face_set_number) const
{
    std::vector<MeshLib::Node*> face_set_nodes;
    std::vector<MeshLib::Element*> face_set_elements;

    for (auto const node : _nodes)
    {
        if (node->isMemberOfFaceSet(face_set_number))
        {
            addFaceSetQuad(node, face_set_number, face_set_nodes,
                           face_set_elements);
        }
    }

    if (face_set_nodes.empty())
    {
        return nullptr;
    }

    for (auto const node : _split_nodes)
    {
        if (node->isMemberOfFaceSet(face_set_number))
        {
            if (node->getAffectedCells()[0])
            {
                addFaceSetQuad(node, face_set_number, face_set_nodes,
                               face_set_elements);
            }
        }
    }

    std::string const mesh_name("FaceSet-" + std::to_string(face_set_number));
    return std::make_unique<MeshLib::Mesh>(mesh_name, face_set_nodes,
                                           face_set_elements);
}

void GocadSGridReader::addFaceSetQuad(
    GocadNode* face_set_node, std::size_t face_set_number,
    std::vector<MeshLib::Node*>& face_set_nodes,
    std::vector<MeshLib::Element*>& face_set_elements) const
{
    std::array<MeshLib::Node*, 4> quad_nodes{};
    quad_nodes[0] = new GocadNode(*face_set_node);
    const std::size_t id(face_set_node->getID());
    std::array<std::size_t, 3> const c(_index_calculator.getCoordsForID(id));

    const FaceDirection dir(face_set_node->getFaceDirection(face_set_number));
    switch (dir)
    {
        case FaceDirection::U:
            quad_nodes[1] = new GocadNode(
                *_nodes[_index_calculator({c[0], c[1] + 1, c[2]})]);
            quad_nodes[2] = new GocadNode(
                *_nodes[_index_calculator({c[0], c[1] + 1, c[2] + 1})]);
            quad_nodes[3] = new GocadNode(
                *_nodes[_index_calculator({c[0], c[1], c[2] + 1})]);
            break;
        case FaceDirection::V:
            quad_nodes[1] = new GocadNode(
                *_nodes[_index_calculator({c[0] + 1, c[1], c[2]})]);
            quad_nodes[2] = new GocadNode(
                *_nodes[_index_calculator({c[0] + 1, c[1], c[2] + 1})]);
            quad_nodes[3] = new GocadNode(
                *_nodes[_index_calculator({c[0], c[1], c[2] + 1})]);
            break;
        case FaceDirection::W:
            quad_nodes[1] = new GocadNode(
                *_nodes[_index_calculator({c[0] + 1, c[1], c[2]})]);
            quad_nodes[2] = new GocadNode(
                *_nodes[_index_calculator({c[0] + 1, c[1] + 1, c[2]})]);
            quad_nodes[3] = new GocadNode(
                *_nodes[_index_calculator({c[0], c[1] + 1, c[2]})]);
            break;
        default:
            ERR("Could not create face for node with id {:d}.", id);
    }
    std::copy(begin(quad_nodes), end(quad_nodes),
              back_inserter(face_set_nodes));
    face_set_elements.push_back(new MeshLib::Quad(quad_nodes));
}

}  // end namespace Gocad
}  // end namespace FileIO
