/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_SECONDARYVARIABLECONTEXT_H
#define PROCESSLIB_SECONDARYVARIABLECONTEXT_H

namespace ProcessLib
{
class SecondaryVariableContext
{
public:
    GlobalIndexType getIndex() const { return _index; }
    void setIndex(GlobalIndexType index) { _index = index; }
private:
    GlobalIndexType _index = 0;
};

}  // namespace ProcessLib

#endif  // PROCESSLIB_SECONDARYVARIABLECONTEXT_H
