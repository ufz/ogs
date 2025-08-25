/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"

// Needs to be exported, see
// https://pybind11.readthedocs.io/en/stable/advanced/misc.html#partitioning-code-over-multiple-extension-modules
class PYBIND11_EXPORT OGSMesh
{
public:
    explicit OGSMesh(MeshLib::Mesh& mesh);

    std::vector<double> getPointCoordinates() const;
    std::pair<std::vector<int>, std::vector<int>> getCells() const;

    std::vector<std::string> getDataArrayNames() const;

    MeshLib::MeshItemType meshItemType(std::string_view const name) const;

    template <typename T>
    pybind11::array& dataArray(std::string const& name)
    {
        auto const data_array_it = data_array_mapping.find(name);
        if (data_array_it != data_array_mapping.end())
        {
            INFO("found data array '{}' with address: {}", name,
                 fmt::ptr(&(data_array_it->second)));
            return data_array_it->second;
        }

        auto& mesh_properties = _mesh.getProperties();

        // check if PropertyVector with specified data type T exists
        auto* property = mesh_properties.getPropertyVector<T>(name);
        if (property == nullptr)
        {
            OGS_FATAL("Couldn't access data array '{}'.", name);
        }

        // Capsule ties lifetime of the numpy array to the PropertyVector
        auto capsule = pybind11::capsule(property, [](void* /*ignored*/) {});

        auto const n_components = property->getNumberOfGlobalComponents();
        if (n_components == 1)
        {
            auto const& [it, success] = data_array_mapping.insert(
                {name,
                 pybind11::array(
                     pybind11::buffer_info(
                         property->data(),  // data pointer
                         sizeof(T),         // item size
                         pybind11::format_descriptor<T>::format(),  // dtype
                                                                    // format
                                                                    // string
                         1,                                         // ndim
                         {property->size()},                        // shape
                         {sizeof(T)}                                // strides
                         ),
                     capsule)});
            if (!success)
            {
                OGS_FATAL("Could not insert data array '{}' into internal map.",
                          name);
            }
            INFO("insert data array '{}' with address: {}", name,
                 fmt::ptr(&(it->second)));
            return it->second;
        }
        else
        {
            // 2D case: (n_items, n_components)
            auto const& [it, success] = data_array_mapping.insert(
                {name,
                 pybind11::array(
                     pybind11::buffer_info(
                         property->data(),  // data pointer
                         sizeof(T),         // item size
                         pybind11::format_descriptor<T>::format(),  // dtype
                                                                    // format
                                                                    // string
                         2,                                         // ndim
                         std::vector<pybind11::ssize_t>{
                             static_cast<pybind11::ssize_t>(property->size() /
                                                            n_components),
                             static_cast<pybind11::ssize_t>(
                                 n_components)},  // shape
                         std::vector<pybind11::ssize_t>{
                             static_cast<pybind11::ssize_t>(sizeof(T) *
                                                            n_components),
                             static_cast<pybind11::ssize_t>(sizeof(T))}
                         // strides
                         ),
                     capsule)});
            if (!success)
            {
                OGS_FATAL("Could not insert data array '{}' into internal map.",
                          name);
            }
            INFO("insert data array '{}' with address: {}", name,
                 fmt::ptr(&(it->second)));
            return it->second;
        }
    }

    pybind11::object dataArray_dispatch(std::string const& name,
                                        std::string const& dtype)
    {
        if (dtype == "double")
        {
            return dataArray<double>(name);
        }
        if (dtype == "float")
        {
            return dataArray<float>(name);
        }
        if (dtype == "int")
        {
            return dataArray<int>(name);
        }
        if (dtype == "int64")
        {
            return dataArray<int64_t>(name);
        }
        if (dtype == "std::size_t")
        {
            return dataArray<std::size_t>(name);
        }
        if (dtype == "char")
        {
            return dataArray<char>(name);
        }

        throw std::runtime_error("Unsupported dtype: " + dtype);
    }

    pybind11::object materialIDs() { return dataArray<int>("MaterialIDs"); }

private:
    MeshLib::Mesh& _mesh;
    std::map<std::string, pybind11::array> data_array_mapping;
};
