/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Implementation of the LisVector class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LisVector.h"

namespace MathLib
{

LisVector::LisVector(std::size_t length)
: _length(length)
{
    lis_vector_create(0, &_vec);
    lis_vector_set_size(_vec, 0, length);
}

LisVector::LisVector(LisVector const &src)
: _length(src.size())
{
	lis_vector_duplicate(src._vec, &_vec);
	lis_vector_copy(src._vec, _vec);
}

LisVector::~LisVector()
{
    lis_vector_destroy(_vec);
}

LisVector& LisVector::operator= (const LisVector &src)
{
	lis_vector_copy(src._vec, _vec);
    return *this;
}

void LisVector::operator+= (const LisVector& v)
{
	lis_vector_axpy(1.0, v._vec, _vec);
}

void LisVector::operator-= (const LisVector& v)
{
	lis_vector_axpy(-1.0, v._vec, _vec);
}

LisVector& LisVector::operator= (double v)
{
    lis_vector_set_all(v, _vec);
    return *this;
}

void LisVector::write (const std::string &filename) const
{
	lis_output_vector(_vec, LIS_FMT_PLAIN, const_cast<char*>(filename.c_str()));
}


} // MathLib

