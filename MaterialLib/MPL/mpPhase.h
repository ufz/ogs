/**
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#ifndef MATERIALLIB_MPL_MPPHASE_H_
#define MATERIALLIB_MPL_MPPHASE_H_

#include "BaseLib/ConfigTree.h"
#include "mpComponent.h"

#include <string>
#include <vector>

namespace MaterialPropertyLib
{
class Phase
{
private:
    std::string _name;
    std::vector<Component*> _components;
    PropertyArray _properties;
public:
    Phase();
    Phase(boost::optional<std::string> const&);
    void createComponents(BaseLib::ConfigTree const&);
    void createProperties(BaseLib::ConfigTree const&);

    Component* component (std::size_t const&);
    std::string name (void);
    std::size_t numberOfComponents(void);

};

} //MaterialPropertyLib



#endif /* MATERIALLIB_MPL_MPPHASE_H_ */
