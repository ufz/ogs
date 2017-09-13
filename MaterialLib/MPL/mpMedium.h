/**
 * \author Norbert Grunwald
 * \date   07.09.2017
 * \brief  
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_MPL_MPMEDIUM_H_
#define MATERIALLIB_MPL_MPMEDIUM_H_

#include "BaseLib/ConfigTree.h"
#include "mpPhase.h"
#include "mpProperty.h"

#include <vector>

namespace MaterialPropertyLib
{
class Medium
{
private:
    std::vector<Phase*> _phases;
    PropertyArray _properties;
public:
    Medium(BaseLib::ConfigTree const&);
    void createPhases (BaseLib::ConfigTree const&);
    void createProperties (BaseLib::ConfigTree const&);

    Phase* phase (std::size_t const);
    std::size_t numberOfPhases(void);

    /// A method that prints out a summary of the constructed medium object
    /// with all its phases, components, and properties. Basically for
    /// debugging purposes only.
    void summary (void);
};

}

#endif /* MATERIALLIB_MPL_MPMEDIUM_H_ */
