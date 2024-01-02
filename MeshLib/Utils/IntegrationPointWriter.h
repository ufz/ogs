/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#include <algorithm>
#include <functional>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#pragma once

namespace MeshLib
{
class Mesh;
class Properties;

struct IntegrationPointWriter final
{
    /// Constructor taking a member function of the \c LocalAssemblerInterface
    /// as a callback.
    ///
    /// The callback \c getIpData must return a \c std::vector<double>.
    template <typename LocalAssemblerInterface, typename... Args>
    IntegrationPointWriter(
        std::string const& name,
        int const n_components,
        int const integration_order,
        std::vector<std::unique_ptr<LocalAssemblerInterface>> const&
            local_assemblers,
        std::vector<double> (LocalAssemblerInterface::*getIpData)(Args...)
            const,
        Args&&... args)
        : _name(name),
          _n_components(n_components),
          _integration_order(integration_order)
    {
        _callback = [&local_assemblers,
                     getIpData,
                     ... f_args = std::forward<Args>(args)]
        {
            // Result containing integration point data for each local
            // assembler.
            std::vector<std::vector<double>> result;
            result.reserve(local_assemblers.size());

            std::transform(begin(local_assemblers), end(local_assemblers),
                           std::back_inserter(result),
                           [&](auto const& la)
                           { return (*la.*getIpData)(f_args...); });

            return result;
        };
    }

    /// Constructor taking an \c accessor function as a callback.
    ///
    /// The accessor function must take a const reference to a
    /// \c LocalAssemblerInterface object as its only argument.
    /// The callback \c accessor must return a \c std::vector<double>.
    template <typename LocalAssemblerInterface, typename Accessor>
    IntegrationPointWriter(
        std::string const& name,
        int const n_components,
        int const integration_order,
        std::vector<std::unique_ptr<LocalAssemblerInterface>> const&
            local_assemblers,
        Accessor accessor)
        : _name(name),
          _n_components(n_components),
          _integration_order(integration_order)
    {
        _callback = [&local_assemblers, accessor]
        {
            // Result containing integration point data for each local
            // assembler.
            std::vector<std::vector<double>> result;
            result.reserve(local_assemblers.size());

            std::transform(begin(local_assemblers), end(local_assemblers),
                           std::back_inserter(result),
                           [&accessor](auto const& la)
                           { return accessor(*la); });

            return result;
        };
    }

    int numberOfComponents() const { return _n_components; }
    int integrationOrder() const { return _integration_order; }
    std::string name() const { return _name; }
    std::vector<std::vector<double>> values() const { return _callback(); }

private:
    std::string const _name;
    int const _n_components;
    int const _integration_order;
    std::function<std::vector<std::vector<double>>()> _callback;
};

/// Add integration point data the the mesh's properties.
///
/// Adds all integration point data arrays given by the input vector and the
/// corresponding meta data as VTK's field data.
/// Integration point data stored as field data (contrary to point or cell
/// data), as plain double arrays. The data is supplemented with information in
/// JSON format, which is stored as array of characters.
void addIntegrationPointDataToMesh(
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
IntegrationPointMetaData getIntegrationPointMetaData(
    MeshLib::Properties const& properties, std::string const& name);
}  // namespace MeshLib
