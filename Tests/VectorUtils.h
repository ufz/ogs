// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <random>
#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "MathLib/LinAlg/UnifiedMatrixSetters.h"

#pragma once

template <typename Vector>
void fillVectorRandomly(Vector& x)
{
    std::random_device rd;
    std::mt19937 random_number_generator(rd());
    std::uniform_real_distribution<double> rnd;

    using Index = typename MathLib::MatrixVectorTraits<Vector>::Index;
    Index const size = x.size();

    for (Index i = 0; i < size; ++i) {
        MathLib::setVector(x, i, rnd(random_number_generator));
    }
#ifdef USE_PETSC
    finalizeVectorAssembly(x);
#endif
}

inline void fillVectorRandomly(std::vector<double>& x)
{
    std::random_device rd;
    std::mt19937 random_number_generator(rd());
    std::uniform_real_distribution<double> rnd;

    for (auto& value : x)
    {
        value = rnd(random_number_generator);
    }
}

template <std::size_t N>
void fillVectorRandomly(std::array<double, N>& x)
{
    std::random_device rd;
    std::mt19937 random_number_generator(rd());
    std::uniform_real_distribution<double> rnd;

    for (auto& value : x)
    {
        value = rnd(random_number_generator);
    }
}
