/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <memory>
#include <nlohmann/json_fwd.hpp>
#include <vector>

#pragma once

namespace MeshLib
{
class Mesh;
}

namespace ProcessLib
{
struct IntegrationPointWriter
{
    virtual ~IntegrationPointWriter() = default;

    virtual int numberOfComponents() const = 0;
    virtual int integrationOrder() const = 0;
    virtual std::string name() const = 0;
    virtual std::vector<std::vector<double>> values() const = 0;
};

/// Add integration point data the the mesh's properties.
///
/// Adds all integration point data arrays given by the input vector and the
/// corresponding meta data as VTK's field data.
/// Integration point data stored as field data (contrary to point or cell
/// data), as plain double arrays. The data is supplemented with information in
/// JSON format, which is stored as array of characters.
void addIntegrationPointWriter(
    MeshLib::Mesh& mesh,
    std::vector<std::unique_ptr<IntegrationPointWriter>> const&
        integration_point_writer);

/// Description of the stored integration point data providing additional
/// information for reconstruction and post-processing.
struct IntegrationPointMetaData
{
    std::string const name;
    int const n_components;
    int const integration_order;
};

/// Returns integration point meta data for the given field name.
///
/// The data is read from a JSON encoded string stored in field data array.
IntegrationPointMetaData getIntegrationPointMetaData(MeshLib::Mesh const& mesh,
                                                     std::string const& name);
}  // namespace ProcessLib
