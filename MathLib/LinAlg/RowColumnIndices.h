// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

namespace MathLib
{

template <typename IDX_TYPE>
struct RowColumnIndices
{
    using LineIndex = typename std::vector<IDX_TYPE>;
    RowColumnIndices(LineIndex const& rows_, LineIndex const& columns_)
        : rows(rows_), columns(columns_)
    {
    }

    LineIndex const& rows;
    LineIndex const& columns;
};

}  // namespace MathLib
