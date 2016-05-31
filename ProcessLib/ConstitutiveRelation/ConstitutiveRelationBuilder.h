/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_CONSTITUTIVERELATION_CONSTITUTIVERELATIONBUILDER_H
#define PROCESSLIB_CONSTITUTIVERELATION_CONSTITUTIVERELATIONBUILDER_H

#include <memory>
#include <type_traits>
#include "ConstitutiveRelation.h"

namespace BaseLib { class ConfigTree; }

namespace ProcessLib {
namespace ConstitutiveRelation {

//! Only used for type erasure.
class ConstitutiveRelationBuilderBase
{
public:
    virtual ~ConstitutiveRelationBuilderBase() = default;
};

template<typename ConstitutiveRelationReturnType,
         typename... ConstitutiveRelationArguments>
class ConstitutiveRelationBuilder
        : public ConstitutiveRelationBuilderBase
{
public:
    virtual std::unique_ptr<
        ConstitutiveRelation<
            ConstitutiveRelationReturnType,
            ConstitutiveRelationArguments...>>
    createConstitutiveRelation(BaseLib::ConfigTree const& config) = 0;
};

} // namespace ConstitutiveRelation
} // namespace ProcessLib

#endif // PROCESSLIB_CONSTITUTIVERELATION_CONSTITUTIVERELATIONBUILDER_H
