/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "SmallDeformationLocalAssemblerFracture.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{
template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
class LocalAssemblerDataFracture final
    : public SmallDeformationLocalAssemblerFracture<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
public:
    LocalAssemblerDataFracture(LocalAssemblerDataFracture const&) = delete;
    LocalAssemblerDataFracture(LocalAssemblerDataFracture&&) = delete;

    LocalAssemblerDataFracture(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        SmallDeformationProcessData<GlobalDim>& process_data)
        : SmallDeformationLocalAssemblerFracture<ShapeFunction,
                                                 IntegrationMethod, GlobalDim>(
              e, local_matrix_size, dofIndex_to_localIndex,
              is_axially_symmetric, integration_order, process_data)
    {
    }
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
