/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_SMALLDEFORMATION_WITH_LIE_SECONDARY_DATA_H_
#define PROCESSLIB_SMALLDEFORMATION_WITH_LIE_SECONDARY_DATA_H_

namespace ProcessLib
{
namespace SmallDeformationWithLIE
{

/// Used by for extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType> N;
};

}  // namespace SmallDeformationWithLIE
}  // namespace ProcessLib

#endif  // PROCESSLIB_SMALLDEFORMATION_WITH_LIE_SECONDARY_DATA_H_
