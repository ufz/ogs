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
#ifndef MATERIALLIB_MPL_MPCOMPONENT_H_
#define MATERIALLIB_MPL_MPCOMPONENT_H_

#include "BaseLib/ConfigTree.h"
#include "mpProperty.h"

namespace MaterialPropertyLib
{
/**
 * \class Component
 * \brief This class defines components (substances).
 * \details The Component class is a base class used for not
 * further specified components. Components are specified by
 * the property 'name'. For specified components we derive
 * special classes from this class (for clarity they are located
 * in the 'components' subfolder).
 */
class Component
{
protected:
    /// The property array of the component.
    PropertyArray _properties;

public:
    Component();
    virtual ~Component() = default;

    /// The method for creating component properties.
    void createProperties(BaseLib::ConfigTree const&);
    /// The method for creating default properties.
    void createDefaultProperties(void);
    /// A get-function for retrieving a cartain property.
    Property* property(PropertyEnum const&);

    void resetPropertyUpdateStatus(void);
};
/*
 * Method for creating a new component based on the specified
 * component name.
 */
std::unique_ptr<Component> newComponent(boost::optional<std::string> const&);

}  // MaterialPropertyLib

#endif /* MATERIALLIB_MPL_MPCOMPONENT_H_ */
