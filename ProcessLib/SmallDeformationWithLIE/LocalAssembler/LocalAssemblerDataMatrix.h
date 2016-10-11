/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_SMALLDEFORMATION_WITH_LIE_LOCALASSEMBLERDATA_H_
#define PROCESSLIB_SMALLDEFORMATION_WITH_LIE_LOCALASSEMBLERDATA_H_

#include "SmallDeformationLocalAssemblerMatrix.h"

namespace ProcessLib
{
namespace SmallDeformationWithLIE
{

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim, int DisplacementDim>
class LocalAssemblerDataMatrix final
    : public SmallDeformationLocalAssemblerMatrix<ShapeFunction, IntegrationMethod,
                                            DisplacementDim>
{
public:
    LocalAssemblerDataMatrix(LocalAssemblerDataMatrix const&) = delete;
    LocalAssemblerDataMatrix(LocalAssemblerDataMatrix&&) = delete;

    LocalAssemblerDataMatrix(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        bool is_axially_symmetric,
        unsigned const integration_order,
        SmallDeformationProcessData<DisplacementDim>& process_data)
        : SmallDeformationLocalAssemblerMatrix<ShapeFunction, IntegrationMethod,
                                         DisplacementDim>(
              e, local_matrix_size, is_axially_symmetric, integration_order,
              process_data)
    {
    }
};

}  // namespace SmallDeformationWithLIE
}  // namespace ProcessLib

#endif  // PROCESSLIB_SMALLDEFORMATION_WITH_LIE_LOCALASSEMBLERDATA_H_
