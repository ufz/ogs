/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_SMALLDEFORMATION_WITH_LIE_LOCALASSEMBLERDATA_FRACTURE_H_
#define PROCESSLIB_SMALLDEFORMATION_WITH_LIE_LOCALASSEMBLERDATA_FRACTURE_H_

#include "SmallDeformationLocalAssemblerFracture.h"

namespace ProcessLib
{
namespace SmallDeformationWithLIE
{

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim, int DisplacementDim>
class LocalAssemblerDataFracture final
    : public SmallDeformationLocalAssemblerFracture<ShapeFunction, IntegrationMethod,
                                            DisplacementDim>
{
public:
    LocalAssemblerDataFracture(LocalAssemblerDataFracture const&) = delete;
    LocalAssemblerDataFracture(LocalAssemblerDataFracture&&) = delete;

    LocalAssemblerDataFracture(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool is_axially_symmetric,
        unsigned const integration_order,
        SmallDeformationProcessData<DisplacementDim>& process_data)
        : SmallDeformationLocalAssemblerFracture<ShapeFunction, IntegrationMethod,
                                         DisplacementDim>(
              e, local_matrix_size, dofIndex_to_localIndex, is_axially_symmetric,
              integration_order, process_data)
    {
    }
};

}  // namespace SmallDeformationWithLIE
}  // namespace ProcessLib

#endif  // PROCESSLIB_SMALLDEFORMATION_WITH_LIE_LOCALASSEMBLERDATA_FRACTURE_H_
