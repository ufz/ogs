/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Implementation of the LisVector class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cassert>

#include "LisVector.h"
#include "LisCheck.h"

namespace MathLib
{
LisVector::LisVector(std::size_t length)
{
    lis_vector_create(0, &_vec);
    lis_vector_set_size(_vec, 0, length);
}

LisVector::LisVector(std::size_t length, double* data)
{
    lis_vector_create(0, &_vec);
    lis_vector_set_size(_vec, 0, length);
    for (std::size_t i = 0; i < length; i++)
        lis_vector_set_value(LIS_INS_VALUE, i, data[i], _vec);
}

LisVector::LisVector(LisVector const& src)
{
    lis_vector_duplicate(src._vec, &_vec);
    lis_vector_copy(src._vec, _vec);
}

LisVector::~LisVector()
{
    lis_vector_destroy(_vec);
}

std::size_t LisVector::size() const
{
    IndexType dummy;
    IndexType size;
    int const ierr = lis_vector_get_size(_vec, &dummy, &size);
    checkLisError(ierr);
    assert(size >= 0);  // For safe implicit conversion to std::size_t.
    return size;
}

void LisVector::write(const std::string& filename) const
{
    lis_output_vector(_vec, LIS_FMT_PLAIN, const_cast<char*>(filename.c_str()));
}

}  // MathLib
