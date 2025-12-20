// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "GetElementRotationMatrices.h"

#include "GetSpaceDimension.h"
#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "MeshLib/Elements/Elements.h"
#include "MeshLib/Mesh.h"

namespace MeshLib
{
std::vector<Eigen::MatrixXd> getElementRotationMatrices(
    int const space_dimension, int const mesh_dimension,
    std::vector<Element*> const& elements)
{
    std::vector<Eigen::MatrixXd> element_rotation_matrices;
    element_rotation_matrices.reserve(elements.size());
    for (auto const* const element : elements)
    {
        int const element_dimension = static_cast<int>(element->getDimension());
        if (element_dimension == space_dimension)
        {
            element_rotation_matrices.emplace_back(Eigen::MatrixXd::Identity(
                element_dimension, element_dimension));
        }
        else
        {
            MeshLib::ElementCoordinatesMappingLocal coordinates_mapping(
                *element, mesh_dimension);

            element_rotation_matrices.emplace_back(
                coordinates_mapping.getRotationMatrixToGlobal().topLeftCorner(
                    space_dimension, element_dimension));
        }
    }

    return element_rotation_matrices;
}
}  // namespace MeshLib
