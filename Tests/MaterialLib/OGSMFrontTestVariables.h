// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
