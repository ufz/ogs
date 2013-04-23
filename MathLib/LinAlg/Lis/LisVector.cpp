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

LisVector::LisVector(unsigned length) : VectorBase<double>(length)
{
    lis_vector_create(0, &_vec);
    lis_vector_set_size(_vec, 0, length);
}

LisVector::~LisVector()
{
    lis_vector_destroy(_vec);
}

void LisVector::printout(std::ostream &os) const
{
    for (unsigned i=0; i<this->getNRows(); i++)
        os << i << " ";
}


} // MathLib

