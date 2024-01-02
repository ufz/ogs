/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

// Currently (09/2022) there are no vectorial and tensorial variables defined
// for use with the extended OGS-MFront interface. But some unit tests need
// them. So we have to provide them here.

#ifdef OGS_USE_MFRONT

#include "MaterialLib/SolidModels/MFront/Variable.h"

struct Vector : MaterialLib::Solids::MFront::Variable<Vector>
{
    constexpr static const char* name = "vector";
    constexpr static mgis::behaviour::Variable::Type type =
        mgis::behaviour::Variable::Type::VECTOR;
};

struct Tensor : MaterialLib::Solids::MFront::Variable<Tensor>
{
    constexpr static const char* name = "tensor";
    constexpr static mgis::behaviour::Variable::Type type =
        mgis::behaviour::Variable::Type::TENSOR;
};
#endif
