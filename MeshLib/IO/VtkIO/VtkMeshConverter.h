// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vtkDataArray.h>
#include <vtkType.h>

#include "BaseLib/DemangleTypeInfo.h"
#include "BaseLib/Logging.h"
#include "GeoLib/Raster.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Properties.h"
#include "MeshLib/PropertyVector.h"

class vtkUnstructuredGrid;  // For conversion vom vtk to ogs mesh
class vtkDataArray;         // For node/cell properties

namespace MeshLib
{
class Mesh;
class Element;
class Node;

/**
 * \brief Converter for VtkUnstructured Grids to OGS meshes
 */
class VtkMeshConverter
{
public:
    /// Converts a vtkUnstructuredGrid object to a Mesh
    static MeshLib::Mesh* convertUnstructuredGrid(
        vtkUnstructuredGrid* grid,
        bool const compute_element_neighbors = false,
        std::string const& mesh_name = "vtkUnstructuredGrid");

private:
    static void convertScalarArrays(vtkUnstructuredGrid& grid,
                                    MeshLib::Mesh& mesh);

    static std::vector<MeshLib::Node*> createNodeVector(
        std::vector<double> const& elevation,
        std::vector<int>& node_idx_map,
        GeoLib::RasterHeader const& header,
        bool use_elevation);

    /// Creates a mesh element vector based on image data
    static std::vector<MeshLib::Element*> createElementVector(
        std::vector<double> const& pix_val,
        std::vector<bool> const& pix_vis,
        std::vector<MeshLib::Node*> const& nodes,
        std::vector<int> const& node_idx_map,
        std::size_t const imgHeight,
        std::size_t const imgWidth,
        MeshElemType elem_type);

    /// Creates a scalar array/mesh property based on pixel values
    template <typename T>
    static void fillPropertyVector(MeshLib::PropertyVector<T>& prop_vec,
                                   std::vector<double> const& pix_val,
                                   std::vector<bool> const& pix_vis,
                                   const std::size_t& imgHeight,
                                   const std::size_t& imgWidth,
                                   MeshElemType elem_type)
    {
        for (std::size_t i = 0; i < imgHeight; i++)
        {
            for (std::size_t j = 0; j < imgWidth; j++)
            {
                if (!pix_vis[i * imgWidth + j])
                {
                    continue;
                }
                T val(static_cast<T>(pix_val[i * (imgWidth + 1) + j]));
                if (elem_type == MeshElemType::TRIANGLE)
                {
                    prop_vec.push_back(val);
                    prop_vec.push_back(val);
                }
                else if (elem_type == MeshElemType::QUAD)
                {
                    prop_vec.push_back(val);
                }
            }
        }
    }

    static void convertArray(vtkDataArray& array,
                             MeshLib::Properties& properties,
                             MeshLib::MeshItemType type);

    template <typename T>
    static void convertTypedArray(vtkDataArray& array,
                                  MeshLib::Properties& properties,
                                  MeshLib::MeshItemType type)
    {
        // Keep for debugging purposes, since the vtk array type and the end
        // type T are not always the same.
        // DBUG(
        //    "Converting a vtk array '{:s}' with data type '{:s}' of "
        //    "size {:d} to a type '{:s}' of size {:d}.",
        //    array.GetName(), array.GetDataTypeAsString(),
        //    array.GetDataTypeSize(), typeid(T).name(), sizeof(T));

        if (sizeof(T) != array.GetDataTypeSize())
        {
            OGS_FATAL(
                "Trying to convert a vtk array '{:s}' with data type '{:s}' of "
                "size {:d} to a different sized type '{:s}' of size {:d}.",
                array.GetName(), array.GetDataTypeAsString(),
                array.GetDataTypeSize(), BaseLib::typeToString<T>(), sizeof(T));
        }

        vtkIdType const nTuples(array.GetNumberOfTuples());
        int const nComponents(array.GetNumberOfComponents());
        char const* const array_name(array.GetName());

        auto* const vec = properties.createNewPropertyVector<T>(
            array_name, type, nTuples, nComponents);
        if (!vec)
        {
            WARN("Array {:s} could not be converted to PropertyVector.",
                 array_name);
            return;
        }
        auto* data_array = static_cast<T*>(array.GetVoidPointer(0));
        vec->assign(std::span<T>(data_array, nTuples * nComponents));
        return;
    }
};

}  // end namespace MeshLib
