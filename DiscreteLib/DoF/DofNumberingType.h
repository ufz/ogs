/**
 * \file   DofNumberingType.h
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief  Helper macros.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#ifndef DOFNUMBERINGTYPE_H_
#define DOFNUMBERINGTYPE_H_

namespace DiscreteLib
{

/**
 * \brief Type of numbering DoF's equation index
 */
struct DofNumberingType
{
    enum type
    {
        BY_VARIABLE,
        BY_POINT
    };
};

} //end

#endif //DOFNUMBERINGTYPE_H_
