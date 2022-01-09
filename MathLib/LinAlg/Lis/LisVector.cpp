/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Implementation of the LisVector class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LisVector.h"

#include "LisCheck.h"

namespace MathLib
{
LisVector::LisVector(std::size_t length)
{
    lis_vector_create(0, &vec_);
    lis_vector_set_size(vec_, 0, length);
}

LisVector::LisVector(std::size_t length, double* data)
{
    lis_vector_create(0, &vec_);
    lis_vector_set_size(vec_, 0, length);
    for (std::size_t i = 0; i < length; i++)
        lis_vector_set_value(LIS_INS_VALUE, i, data[i], vec_);
}

LisVector::LisVector(LisVector const& src)
{
    lis_vector_duplicate(src.vec_, &vec_);
    lis_vector_copy(src.vec_, vec_);
}

LisVector::~LisVector()
{
    lis_vector_destroy(vec_);
}

std::size_t LisVector::size() const
{
    IndexType dummy;
    IndexType size;
    int const ierr = lis_vector_get_size(vec_, &dummy, &size);
    checkLisError(ierr);
    assert(size >= 0);  // For safe implicit conversion to std::size_t.
    return size;
}

void LisVector::write(const std::string& filename) const
{
    lis_output_vector(vec_, LIS_FMT_PLAIN, const_cast<char*>(filename.c_str()));
}

}  // namespace MathLib
