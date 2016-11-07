/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_LIE_SMALLDEFORMATION_LOCALASSEMBLERDATAMATRIXNEARFRACTURE_H_
#define PROCESSLIB_LIE_SMALLDEFORMATION_LOCALASSEMBLERDATAMATRIXNEARFRACTURE_H_

#include "SmallDeformationLocalAssemblerMatrixNearFracture.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim, int DisplacementDim>
class LocalAssemblerDataMatrixNearFracture final
    : public SmallDeformationLocalAssemblerMatrixNearFracture<ShapeFunction, IntegrationMethod,
                                            DisplacementDim>
{
public:
    LocalAssemblerDataMatrixNearFracture(LocalAssemblerDataMatrixNearFracture const&) = delete;
    LocalAssemblerDataMatrixNearFracture(LocalAssemblerDataMatrixNearFracture&&) = delete;

    LocalAssemblerDataMatrixNearFracture(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        SmallDeformationProcessData<DisplacementDim>& process_data)
        : SmallDeformationLocalAssemblerMatrixNearFracture<ShapeFunction, IntegrationMethod,
                                         DisplacementDim>(
              e, n_variables, local_matrix_size, dofIndex_to_localIndex,
              is_axially_symmetric, integration_order, process_data)
    {
    }
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib

#endif // PROCESSLIB_LIE_SMALLDEFORMATION_LOCALASSEMBLERDATAMATRIXNEARFRACTURE_H_
