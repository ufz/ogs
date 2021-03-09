/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <spdlog/spdlog.h>

#include <boost/dynamic_bitset.hpp>
#include <cstddef>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "CoordinateSystem.h"
#include "GocadNode.h"
#include "IndexCalculator.h"
#include "Layer.h"
#include "MeshLib/Elements/Element.h"
#include "Property.h"
#include "Region.h"

namespace FileIO
{
namespace Gocad
{
class GocadSGridReader final
{
public:
    /**
     * Constructor takes as argument the Gocad .sg text file.
     * @param fname file name
     */
    explicit GocadSGridReader(std::string const& fname);
    ~GocadSGridReader();

    GocadSGridReader() = delete;
    GocadSGridReader(GocadSGridReader&& src) = delete;
    GocadSGridReader(GocadSGridReader const& src) = delete;
    GocadSGridReader& operator=(GocadSGridReader&& rhs) = delete;
    GocadSGridReader& operator=(GocadSGridReader const& rhs) = delete;

    std::unique_ptr<MeshLib::Mesh> getMesh() const;
    std::unique_ptr<MeshLib::Mesh> getFaceSetMesh(
        std::size_t const face_set_number) const;

    std::vector<std::string> getPropertyNames() const;

private:
    using Bitset = boost::dynamic_bitset<>;

    void parseDims(std::string const& line);
    void parseFileName(std::string const& line,
                       std::string& result_string) const;
    void parseHeader(std::istream& in);
    void parseFaceSet(std::string& line, std::istream& in);

    void readNodesBinary();
    std::vector<Bitset> readRegionFlagsBinary() const;
    void readElementPropertiesBinary();
    void mapRegionFlagsToCellProperties(std::vector<Bitset> const& rf);

    std::vector<MeshLib::Element*> createElements(
        std::vector<MeshLib::Node*> const& nodes) const;

    // split handling
    void readSplitInformation();
    void applySplitInformation(
        std::vector<MeshLib::Node*>& nodes,
        std::vector<MeshLib::Element*> const& elements) const;
    static void modifyElement(MeshLib::Element* hex,
                              MeshLib::Node const* node2sub,
                              MeshLib::Node* substitute_node);

    void addFaceSetQuad(
        GocadNode* face_set_node, std::size_t face_set_number,
        std::vector<MeshLib::Node*>& face_set_nodes,
        std::vector<MeshLib::Element*>& face_set_elements) const;

    std::optional<Gocad::Property const&> getProperty(
        std::string const& name) const;
    void addGocadPropertiesToMesh(MeshLib::Mesh& mesh) const;

    std::string const& _fname;
    std::string const _path;
    // data read from sg file
    Gocad::IndexCalculator _index_calculator;
    Gocad::CoordinateSystem _coordinate_system;
    std::string _pnts_fname;
    std::string _flags_fname;
    std::string _region_flags_fname;

    std::vector<Gocad::Region> regions;
    std::vector<Gocad::Layer> layers;
    std::size_t _n_face_sets;

    bool _double_precision_binary;
    bool _bin_pnts_in_double_precision;

    // data read from binary points file
    std::vector<GocadNode*> _nodes;
    std::vector<GocadSplitNode*> _split_nodes;
    // properties
    std::vector<Gocad::Property> _property_meta_data_vecs;
};  // end class GocadSGridReader

}  // end namespace Gocad
}  // end namespace FileIO
