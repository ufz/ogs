/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#include <functional>
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
struct IntegrationPointWriter final
{
    IntegrationPointWriter(std::string const& name,
                           int const n_components,
                           int const integration_order,
                           std::function<std::vector<std::vector<double>>()>
                               callback)
        : name_(name),
          n_components_(n_components),
          integration_order_(integration_order),
          callback_(callback)
    {
    }

    int numberOfComponents() const { return n_components_; }
    int integrationOrder() const { return integration_order_; }
    std::string name() const { return name_; }
    std::vector<std::vector<double>> values() const { return callback_(); }

private:
    std::string const name_;
    int const n_components_;
    int const integration_order_;
    std::function<std::vector<std::vector<double>>()> callback_;
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
