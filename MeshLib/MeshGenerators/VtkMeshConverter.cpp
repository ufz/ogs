/**
 * \file
 * \author Karsten Rink
 * \date   2011-08-23
 * \brief  Implementation of the VtkMeshConverter class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkMeshConverter.h"

#include <cstdint>

#include "MeshLib/Elements/Elements.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

// Conversion from Image to QuadMesh
#include <vtkBitArray.h>
#include <vtkCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkLongArray.h>
#include <vtkLongLongArray.h>
#include <vtkPointData.h>
#include <vtkShortArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedLongArray.h>
#include <vtkUnsignedLongLongArray.h>
#include <vtkUnsignedShortArray.h>

// Conversion from vtkUnstructuredGrid
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkUnsignedIntArray.h>
#include <vtkUnstructuredGrid.h>

namespace MeshLib
{
namespace detail
{
template <class T_ELEMENT>
MeshLib::Element* createElementWithSameNodeOrder(
    const std::vector<MeshLib::Node*>& nodes, vtkIdList* const node_ids,
    std::size_t const element_id)
{
    auto** ele_nodes = new MeshLib::Node*[T_ELEMENT::n_all_nodes];
    for (unsigned k(0); k < T_ELEMENT::n_all_nodes; k++)
    {
        ele_nodes[k] = nodes[node_ids->GetId(k)];
    }
    return new T_ELEMENT(ele_nodes, element_id);
}
}  // namespace detail

MeshLib::Mesh* VtkMeshConverter::convertUnstructuredGrid(
    vtkUnstructuredGrid* grid, std::string const& mesh_name)
{
    if (!grid)
    {
        return nullptr;
    }

    // set mesh nodes
    const std::size_t nNodes = grid->GetPoints()->GetNumberOfPoints();
    std::vector<MeshLib::Node*> nodes(nNodes);
    double* coords = nullptr;
    for (std::size_t i = 0; i < nNodes; i++)
    {
        coords = grid->GetPoints()->GetPoint(i);
        nodes[i] = new MeshLib::Node(coords[0], coords[1], coords[2], i);
    }

    // set mesh elements
    const std::size_t nElems = grid->GetNumberOfCells();
    std::vector<MeshLib::Element*> elements(nElems);
    auto node_ids = vtkSmartPointer<vtkIdList>::New();
    for (std::size_t i = 0; i < nElems; i++)
    {
        MeshLib::Element* elem;
        grid->GetCellPoints(i, node_ids);

        int cell_type = grid->GetCellType(i);
        switch (cell_type)
        {
            case VTK_VERTEX:
            {
                elem = detail::createElementWithSameNodeOrder<MeshLib::Point>(
                    nodes, node_ids, i);
                break;
            }
            case VTK_LINE:
            {
                elem = detail::createElementWithSameNodeOrder<MeshLib::Line>(
                    nodes, node_ids, i);
                break;
            }
            case VTK_TRIANGLE:
            {
                elem = detail::createElementWithSameNodeOrder<MeshLib::Tri>(
                    nodes, node_ids, i);
                break;
            }
            case VTK_QUAD:
            {
                elem = detail::createElementWithSameNodeOrder<MeshLib::Quad>(
                    nodes, node_ids, i);
                break;
            }
            case VTK_PIXEL:
            {
                auto** quad_nodes = new MeshLib::Node*[4];
                quad_nodes[0] = nodes[node_ids->GetId(0)];
                quad_nodes[1] = nodes[node_ids->GetId(1)];
                quad_nodes[2] = nodes[node_ids->GetId(3)];
                quad_nodes[3] = nodes[node_ids->GetId(2)];
                elem = new MeshLib::Quad(quad_nodes, i);
                break;
            }
            case VTK_TETRA:
            {
                elem = detail::createElementWithSameNodeOrder<MeshLib::Tet>(
                    nodes, node_ids, i);
                break;
            }
            case VTK_HEXAHEDRON:
            {
                elem = detail::createElementWithSameNodeOrder<MeshLib::Hex>(
                    nodes, node_ids, i);
                break;
            }
            case VTK_VOXEL:
            {
                auto** voxel_nodes = new MeshLib::Node*[8];
                voxel_nodes[0] = nodes[node_ids->GetId(0)];
                voxel_nodes[1] = nodes[node_ids->GetId(1)];
                voxel_nodes[2] = nodes[node_ids->GetId(3)];
                voxel_nodes[3] = nodes[node_ids->GetId(2)];
                voxel_nodes[4] = nodes[node_ids->GetId(4)];
                voxel_nodes[5] = nodes[node_ids->GetId(5)];
                voxel_nodes[6] = nodes[node_ids->GetId(7)];
                voxel_nodes[7] = nodes[node_ids->GetId(6)];
                elem = new MeshLib::Hex(voxel_nodes, i);
                break;
            }
            case VTK_PYRAMID:
            {
                elem = detail::createElementWithSameNodeOrder<MeshLib::Pyramid>(
                    nodes, node_ids, i);
                break;
            }
            case VTK_WEDGE:
            {
                auto** prism_nodes = new MeshLib::Node*[6];
                for (unsigned j = 0; j < 3; ++j)
                {
                    prism_nodes[j] = nodes[node_ids->GetId(j + 3)];
                    prism_nodes[j + 3] = nodes[node_ids->GetId(j)];
                }
                elem = new MeshLib::Prism(prism_nodes, i);
                break;
            }
            case VTK_QUADRATIC_EDGE:
            {
                elem = detail::createElementWithSameNodeOrder<MeshLib::Line3>(
                    nodes, node_ids, i);
                break;
            }
            case VTK_QUADRATIC_TRIANGLE:
            {
                elem = detail::createElementWithSameNodeOrder<MeshLib::Tri6>(
                    nodes, node_ids, i);
                break;
            }
            case VTK_QUADRATIC_QUAD:
            {
                elem = detail::createElementWithSameNodeOrder<MeshLib::Quad8>(
                    nodes, node_ids, i);
                break;
            }
            case VTK_BIQUADRATIC_QUAD:
            {
                elem = detail::createElementWithSameNodeOrder<MeshLib::Quad9>(
                    nodes, node_ids, i);
                break;
            }
            case VTK_QUADRATIC_TETRA:
            {
                elem = detail::createElementWithSameNodeOrder<MeshLib::Tet10>(
                    nodes, node_ids, i);
                break;
            }
            case VTK_QUADRATIC_HEXAHEDRON:
            {
                elem = detail::createElementWithSameNodeOrder<MeshLib::Hex20>(
                    nodes, node_ids, i);
                break;
            }
            case VTK_QUADRATIC_PYRAMID:
            {
                elem =
                    detail::createElementWithSameNodeOrder<MeshLib::Pyramid13>(
                        nodes, node_ids, i);
                break;
            }
            case VTK_QUADRATIC_WEDGE:
            {
                elem = detail::createElementWithSameNodeOrder<MeshLib::Prism15>(
                    nodes, node_ids, i);

                break;
            }
            default:
                ERR("VtkMeshConverter::convertUnstructuredGrid(): Unknown mesh "
                    "element type '{:d}'.",
                    cell_type);
                return nullptr;
        }

        elements[i] = elem;
    }

    MeshLib::Mesh* mesh = new MeshLib::Mesh(mesh_name, nodes, elements);
    convertScalarArrays(*grid, *mesh);

    return mesh;
}

void VtkMeshConverter::convertScalarArrays(vtkUnstructuredGrid& grid,
                                           MeshLib::Mesh& mesh)
{
    vtkPointData* point_data = grid.GetPointData();
    auto const n_point_arrays =
        static_cast<unsigned>(point_data->GetNumberOfArrays());
    for (unsigned i = 0; i < n_point_arrays; ++i)
    {
        convertArray(*point_data->GetArray(i),
                     mesh.getProperties(),
                     MeshLib::MeshItemType::Node);
    }

    vtkCellData* cell_data = grid.GetCellData();
    auto const n_cell_arrays =
        static_cast<unsigned>(cell_data->GetNumberOfArrays());
    for (unsigned i = 0; i < n_cell_arrays; ++i)
    {
        convertArray(*cell_data->GetArray(i),
                     mesh.getProperties(),
                     MeshLib::MeshItemType::Cell);
    }

    vtkFieldData* field_data = grid.GetFieldData();
    auto const n_field_arrays =
        static_cast<unsigned>(field_data->GetNumberOfArrays());
    for (unsigned i = 0; i < n_field_arrays; ++i)
    {
        convertArray(
            *vtkDataArray::SafeDownCast(field_data->GetAbstractArray(i)),
            mesh.getProperties(),
            MeshLib::MeshItemType::IntegrationPoint);
    }
}

void VtkMeshConverter::convertArray(vtkDataArray& array,
                                    MeshLib::Properties& properties,
                                    MeshLib::MeshItemType type)
{
    if (vtkDoubleArray::SafeDownCast(&array))
    {
        VtkMeshConverter::convertTypedArray<double>(array, properties, type);
        return;
    }

    if (vtkFloatArray::SafeDownCast(&array))
    {
        VtkMeshConverter::convertTypedArray<float>(array, properties, type);
        return;
    }

    if (vtkBitArray::SafeDownCast(&array))
    {
        VtkMeshConverter::convertTypedArray<bool>(array, properties, type);
        return;
    }

    // This also covers the vtkTypeInt8Array type, which is derived from the
    // vtkCharArray type.
    if (vtkCharArray::SafeDownCast(&array))
    {
        if constexpr (sizeof(std::int8_t) != sizeof(char))
        {
            OGS_FATAL(
                "Array '{:s}' in VTU file uses unsupported data type '{:s}' "
                "not convertible to char.",
                array.GetName(), array.GetDataTypeAsString());
        }
        VtkMeshConverter::convertTypedArray<char>(array, properties, type);
        return;
    }

    // This also covers the vtkTypeInt16Array type, which is derived from the
    // vtkShortArray type.
    if (vtkShortArray::SafeDownCast(&array))
    {
        if constexpr (sizeof(std::int16_t) != sizeof(short))
        {
            OGS_FATAL(
                "Array '{:s}' in VTU file uses unsupported data type '{:s}' "
                "not convertible to short.",
                array.GetName(), array.GetDataTypeAsString());
        }
        VtkMeshConverter::convertTypedArray<short>(array, properties, type);
        return;
    }

    // This also covers the vtkTypeInt32Array type, which is derived from the
    // vtkIntArray type.
    if (vtkIntArray::SafeDownCast(&array))
    {
        if constexpr (sizeof(std::int32_t) != sizeof(int))
        {
            OGS_FATAL(
                "Array '{:s}' in VTU file uses unsupported data type '{:s}' "
                "not convertible to int.",
                array.GetName(), array.GetDataTypeAsString());
        }
        VtkMeshConverter::convertTypedArray<int>(array, properties, type);
        return;
    }

    // This is not a unique conversion depending on the sizes of int, long, and
    // long long. Converting to the apparently smaller type.
    if (vtkLongArray::SafeDownCast(&array))
    {
        if constexpr (sizeof(long) == sizeof(int))
        {
            VtkMeshConverter::convertTypedArray<int>(array, properties, type);
            return;
        }
        if constexpr (sizeof(long long) == sizeof(long))
        {
            VtkMeshConverter::convertTypedArray<long>(array, properties, type);
            return;
        }

        OGS_FATAL(
            "Array '{:s}' in VTU file uses unsupported data type '{:s}' ({:d}) "
            "not convertible to int ({:d}), or long ({:d}) types.",
            array.GetName(), array.GetDataTypeAsString(),
            array.GetDataTypeSize(), sizeof(int), sizeof(long));
    }

    // Ensure, vtkTypeInt64Array and vtkLongLongArray are using the same type.
    static_assert(sizeof(std::int64_t) == sizeof(long long));

    // This also covers the vtkTypeInt64Array type, which is derived from the
    // vtkLongLongArray type.
    // Converting to the apparently smaller type.
    if (vtkLongLongArray::SafeDownCast(&array))
    {
        if constexpr (sizeof(long long) == sizeof(long))
        {
            VtkMeshConverter::convertTypedArray<long>(array, properties, type);
            return;
        }
        else  // ll > l. Other cases are not possible because ll >= l in c++.
        {
            VtkMeshConverter::convertTypedArray<long long>(array, properties,
                                                           type);
            return;
        }
    }

    // This also covers the vtkTypeUInt8Array type, which is derived from the
    // vtkUnsignedCharArray type.
    if (vtkUnsignedCharArray::SafeDownCast(&array))
    {
        if constexpr (sizeof(std::uint8_t) != sizeof(unsigned char))
        {
            OGS_FATAL(
                "Array '{:s}' in VTU file uses unsupported data type '{:s}' "
                "not convertible to unsigned char.",
                array.GetName(), array.GetDataTypeAsString());
        }
        VtkMeshConverter::convertTypedArray<unsigned char>(array, properties,
                                                           type);
        return;
    }

    // This also covers the vtkTypeUInt16Array type, which is derived from the
    // vtkUnsignedShortArray type.
    if (vtkUnsignedShortArray::SafeDownCast(&array))
    {
        if constexpr (sizeof(std::uint16_t) != sizeof(unsigned short))
        {
            OGS_FATAL(
                "Array '{:s}' in VTU file uses unsupported data type '{:s}' "
                "not convertible to unsigned short.",
                array.GetName(), array.GetDataTypeAsString());
        }
        VtkMeshConverter::convertTypedArray<unsigned short>(array, properties,
                                                            type);
        return;
    }

    // This also covers the vtkTypeUInt32Array type, which is derived from the
    // vtkUnsignedIntArray type.
    if (vtkUnsignedIntArray::SafeDownCast(&array))
    {
        if constexpr (sizeof(std::uint32_t) != sizeof(unsigned))
        {
            OGS_FATAL(
                "Array '{:s}' in VTU file uses unsupported data type '{:s}' "
                "not convertible to unsigned.",
                array.GetName(), array.GetDataTypeAsString());
        }

        // MaterialIDs are assumed to be integers
        if (std::strncmp(array.GetName(), "MaterialIDs", 11) == 0)
        {
            VtkMeshConverter::convertTypedArray<int>(array, properties, type);
        }
        else
        {
            VtkMeshConverter::convertTypedArray<unsigned>(array, properties,
                                                          type);
        }

        return;
    }

    // This is not a unique conversion depending on the sizes of unsigned,
    // unsigned long, and unsigned long long. Converting to the apparently
    // smaller type.
    if (vtkUnsignedLongArray::SafeDownCast(&array))
    {
        if constexpr (sizeof(unsigned long) == sizeof(unsigned))
        {
            VtkMeshConverter::convertTypedArray<unsigned>(array, properties,
                                                          type);
            return;
        }
        if constexpr (sizeof(unsigned long long) == sizeof(unsigned long))
        {
            VtkMeshConverter::convertTypedArray<unsigned long>(
                array, properties, type);
            return;
        }

        OGS_FATAL(
            "Array '{:s}' in VTU file uses unsupported data type '{:s}' ({:d}) "
            "not convertible to unsigned ({:d}), or unsigned long ({:d}) "
            "types.",
            array.GetName(), array.GetDataTypeAsString(),
            array.GetDataTypeSize(), sizeof(unsigned), sizeof(unsigned long));
    }

    // Ensure, vtkTypeUInt64Array and vtkUnsignedLongLongArray are using the
    // same type.
    static_assert(sizeof(std::uint64_t) == sizeof(unsigned long long));

    // This also covers the vtkTypeUInt64Array type, which is derived from the
    // vtkUnsignedLongLongArray type.
    // Converting to the apparently smaller type.
    if (vtkUnsignedLongLongArray::SafeDownCast(&array))
    {
        if constexpr (sizeof(unsigned long long) == sizeof(unsigned long))
        {
            VtkMeshConverter::convertTypedArray<unsigned long>(
                array, properties, type);
            return;
        }
        else  // ull > ul. Other cases are not possible because ull >= ul in
              // c++.
        {
            VtkMeshConverter::convertTypedArray<unsigned long long>(
                array, properties, type);
            return;
        }
    }

    WARN(
        "Array '{:s}' in VTU file uses unsupported data type '{:s}' of size "
        "{:d}. The data array will not be available.",
        array.GetName(), array.GetDataTypeAsString(), array.GetDataTypeSize());
}

}  // end namespace MeshLib
