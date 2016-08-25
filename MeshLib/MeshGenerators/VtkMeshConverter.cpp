/**
 * \file
 * \author Karsten Rink
 * \date   2011-08-23
 * \brief  Implementation of the VtkMeshConverter class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkMeshConverter.h"

#include "MeshLib/Elements/Elements.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

// Conversion from Image to QuadMesh
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkBitArray.h>
#include <vtkCharArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedLongArray.h>
#include <vtkUnsignedLongLongArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

// Conversion from vtkUnstructuredGrid
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnsignedIntArray.h>

namespace MeshLib
{

namespace detail
{
template <class T_ELEMENT>
MeshLib::Element* createElementWithSameNodeOrder(const std::vector<MeshLib::Node*> &nodes,
        vtkIdList* const node_ids)
{
    MeshLib::Node** ele_nodes = new MeshLib::Node*[T_ELEMENT::n_all_nodes];
    for (unsigned k(0); k<T_ELEMENT::n_all_nodes; k++)
        ele_nodes[k] = nodes[node_ids->GetId(k)];
    return new T_ELEMENT(ele_nodes);
}
}

MeshLib::Mesh* VtkMeshConverter::convertUnstructuredGrid(vtkUnstructuredGrid* grid, std::string const& mesh_name)
{
    if (!grid)
        return nullptr;

    // set mesh nodes
    const std::size_t nNodes = grid->GetPoints()->GetNumberOfPoints();
    std::vector<MeshLib::Node*> nodes(nNodes);
    double* coords = nullptr;
    for (std::size_t i = 0; i < nNodes; i++)
    {
        coords = grid->GetPoints()->GetPoint(i);
        nodes[i] = new MeshLib::Node(coords[0], coords[1], coords[2]);
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
        case VTK_LINE: {
            elem = detail::createElementWithSameNodeOrder<MeshLib::Line>(nodes, node_ids);
            break;
        }
        case VTK_TRIANGLE: {
            elem = detail::createElementWithSameNodeOrder<MeshLib::Tri>(nodes, node_ids);
            break;
        }
        case VTK_QUAD: {
            elem = detail::createElementWithSameNodeOrder<MeshLib::Quad>(nodes, node_ids);
            break;
        }
        case VTK_PIXEL: {
            MeshLib::Node** quad_nodes = new MeshLib::Node*[4];
            quad_nodes[0] = nodes[node_ids->GetId(0)];
            quad_nodes[1] = nodes[node_ids->GetId(1)];
            quad_nodes[2] = nodes[node_ids->GetId(3)];
            quad_nodes[3] = nodes[node_ids->GetId(2)];
            elem = new MeshLib::Quad(quad_nodes);
            break;
        }
        case VTK_TETRA: {
            elem = detail::createElementWithSameNodeOrder<MeshLib::Tet>(nodes, node_ids);
            break;
        }
        case VTK_HEXAHEDRON: {
            elem = detail::createElementWithSameNodeOrder<MeshLib::Hex>(nodes, node_ids);
            break;
        }
        case VTK_VOXEL: {
            MeshLib::Node** voxel_nodes = new MeshLib::Node*[8];
            voxel_nodes[0] = nodes[node_ids->GetId(0)];
            voxel_nodes[1] = nodes[node_ids->GetId(1)];
            voxel_nodes[2] = nodes[node_ids->GetId(3)];
            voxel_nodes[3] = nodes[node_ids->GetId(2)];
            voxel_nodes[4] = nodes[node_ids->GetId(4)];
            voxel_nodes[5] = nodes[node_ids->GetId(5)];
            voxel_nodes[6] = nodes[node_ids->GetId(7)];
            voxel_nodes[7] = nodes[node_ids->GetId(6)];
            elem = new MeshLib::Hex(voxel_nodes);
            break;
        }
        case VTK_PYRAMID: {
            elem = detail::createElementWithSameNodeOrder<MeshLib::Pyramid>(nodes, node_ids);
            break;
        }
        case VTK_WEDGE: {
            MeshLib::Node** prism_nodes = new MeshLib::Node*[6];
            for (unsigned i=0; i<3; ++i)
            {
                prism_nodes[i] = nodes[node_ids->GetId(i+3)];
                prism_nodes[i+3] = nodes[node_ids->GetId(i)];
            }
            elem = new MeshLib::Prism(prism_nodes);
            break;
        }
        case VTK_QUADRATIC_EDGE: {
            elem = detail::createElementWithSameNodeOrder<MeshLib::Line3>(nodes, node_ids);
            break;
        }
        case VTK_QUADRATIC_TRIANGLE: {
            elem = detail::createElementWithSameNodeOrder<MeshLib::Tri6>(nodes, node_ids);
            break;
        }
        case VTK_QUADRATIC_QUAD: {
            elem = detail::createElementWithSameNodeOrder<MeshLib::Quad8>(nodes, node_ids);
            break;
        }
        case VTK_BIQUADRATIC_QUAD: {
            elem = detail::createElementWithSameNodeOrder<MeshLib::Quad9>(nodes, node_ids);
            break;
        }
        case VTK_QUADRATIC_TETRA: {
            elem = detail::createElementWithSameNodeOrder<MeshLib::Tet10>(nodes, node_ids);
            break;
        }
        case VTK_QUADRATIC_HEXAHEDRON: {
            elem = detail::createElementWithSameNodeOrder<MeshLib::Hex20>(nodes, node_ids);
            break;
        }
        case VTK_QUADRATIC_PYRAMID: {
            elem = detail::createElementWithSameNodeOrder<MeshLib::Pyramid13>(nodes, node_ids);
            break;
        }
        case VTK_QUADRATIC_WEDGE: {
            MeshLib::Node** prism_nodes = new MeshLib::Node*[15];
            for (unsigned i=0; i<3; ++i)
            {
                prism_nodes[i] = nodes[node_ids->GetId(i+3)];
                prism_nodes[i+3] = nodes[node_ids->GetId(i)];
            }
            for (unsigned i=0; i<3; ++i)
                prism_nodes[6+i] = nodes[node_ids->GetId(8-i)];
            prism_nodes[9] = nodes[node_ids->GetId(12)];
            prism_nodes[10] = nodes[node_ids->GetId(14)];
            prism_nodes[11] = nodes[node_ids->GetId(13)];
            for (unsigned i=0; i<3; ++i)
                prism_nodes[12+i] = nodes[node_ids->GetId(11-i)];
            elem = new MeshLib::Prism15(prism_nodes);
            break;
        }
        default:
            ERR("VtkMeshConverter::convertUnstructuredGrid(): Unknown mesh element type \"%d\".", cell_type);
            return nullptr;
        }

        elements[i] = elem;
    }

    MeshLib::Mesh* mesh = new MeshLib::Mesh(mesh_name, nodes, elements);
    convertScalarArrays(*grid, *mesh);

    return mesh;
}

void VtkMeshConverter::convertScalarArrays(vtkUnstructuredGrid &grid, MeshLib::Mesh &mesh)
{
    vtkPointData* point_data = grid.GetPointData();
    unsigned const n_point_arrays = static_cast<unsigned>(point_data->GetNumberOfArrays());
    for (unsigned i=0; i<n_point_arrays; ++i)
        convertArray(*point_data->GetArray(i), mesh.getProperties(), MeshLib::MeshItemType::Node);

    vtkCellData* cell_data = grid.GetCellData();
    unsigned const n_cell_arrays = static_cast<unsigned>(cell_data->GetNumberOfArrays());
    for (unsigned i=0; i<n_cell_arrays; ++i)
        convertArray(*cell_data->GetArray(i), mesh.getProperties(), MeshLib::MeshItemType::Cell);
}

void VtkMeshConverter::convertArray(vtkDataArray &array, MeshLib::Properties &properties, MeshLib::MeshItemType type)
{
    if (vtkDoubleArray::SafeDownCast(&array))
    {
        VtkMeshConverter::convertTypedArray<double>(array, properties, type);
        return;
    }

    if (vtkIntArray::SafeDownCast(&array))
    {
        VtkMeshConverter::convertTypedArray<int>(array, properties, type);
        return;
    }

    if (vtkBitArray::SafeDownCast(&array))
    {
        VtkMeshConverter::convertTypedArray<bool>(array, properties, type);
        return;
    }

    if (vtkCharArray::SafeDownCast(&array))
    {
        VtkMeshConverter::convertTypedArray<char>(array, properties, type);
        return;
    }

    if (vtkUnsignedLongArray::SafeDownCast(&array))
    {
        VtkMeshConverter::convertTypedArray<unsigned long>(array, properties, type);
        return;
    }

    if (vtkUnsignedLongLongArray::SafeDownCast(&array))
    {
        VtkMeshConverter::convertTypedArray<unsigned long long>(
            array, properties, type);
        return;
    }

    if (vtkUnsignedIntArray::SafeDownCast(&array))
    {
        // MaterialIDs are assumed to be integers
        if(std::strncmp(array.GetName(), "MaterialIDs", 11) == 0)
            VtkMeshConverter::convertTypedArray<int>(array, properties, type);
        else
            VtkMeshConverter::convertTypedArray<unsigned>(array, properties, type);

        return;
    }

    ERR ("Array \"%s\" in VTU file uses unsupported data type.", array.GetName());
    return;
}

} // end namespace MeshLib
