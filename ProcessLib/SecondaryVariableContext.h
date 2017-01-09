/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_SECONDARYVARIABLECONTEXT_H
#define PROCESSLIB_SECONDARYVARIABLECONTEXT_H

namespace ProcessLib
{
/*! Provides a "global" context, e.g., for NamedFunction's.
 *
 * Currently the context comprises an index. It can be used, e.g. to compute a
 * NamedFunction which needs additional external data apart from the unbound
 * variables that are configured.
 */
struct SecondaryVariableContext
{
public:
    //! Points to the position in a GlobalVector being read or written right
    //! now.
    GlobalIndexType index = 0;
};

}  // namespace ProcessLib

#endif  // PROCESSLIB_SECONDARYVARIABLECONTEXT_H
