/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "SmallDeformationLocalAssemblerMatrix.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{
template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
class LocalAssemblerDataMatrix final
    : public SmallDeformationLocalAssemblerMatrix<ShapeFunction,
                                                  IntegrationMethod, GlobalDim>
{
public:
    LocalAssemblerDataMatrix(LocalAssemblerDataMatrix const&) = delete;
    LocalAssemblerDataMatrix(LocalAssemblerDataMatrix&&) = delete;

    LocalAssemblerDataMatrix(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        SmallDeformationProcessData<GlobalDim>& process_data)
        : SmallDeformationLocalAssemblerMatrix<ShapeFunction, IntegrationMethod,
                                               GlobalDim>(
              e, local_matrix_size, is_axially_symmetric, integration_order,
              process_data)
    {
    }
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
