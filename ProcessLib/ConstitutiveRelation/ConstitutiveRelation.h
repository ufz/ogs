/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_CONSTITUTIVERELATION_CONSTITUTIVERELATION_H
#define PROCESSLIB_CONSTITUTIVERELATION_CONSTITUTIVERELATION_H

#include <type_traits>

#include "ProcessLib/Parameter.h" // TODO remove

namespace ProcessLib {
namespace ConstitutiveRelation {

template<typename ReturnType, typename... Arguments>
class ConstitutiveRelation
{
public:
    virtual ReturnType getValue(
        double const t,
        double const*const x,
        GlobalIndexType const node,
        MeshLib::Element const& element,
        std::size_t const integration_point,
        Arguments const&... arguments) const = 0;

    // TODO Default implementation? Central differences?
    virtual
    typename std::enable_if<sizeof...(Arguments)!=0, ReturnType>::type
    getDerivative(
        unsigned const derivative_in_direction_of_argument,
        double const t,
        double const*const x,
        GlobalIndexType const node,
        MeshLib::Element const& element,
        std::size_t const integration_point,
        Arguments const&... arguments) const = 0;

    virtual ~ConstitutiveRelation() = default;
};

} // namespace ConstitutiveRelation
} // namespace ProcessLib

#endif // PROCESSLIB_CONSTITUTIVERELATION_CONSTITUTIVERELATION_H
