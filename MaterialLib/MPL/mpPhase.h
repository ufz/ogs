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
#include "mpProperty.h"

#include <string>
#include <vector>

namespace MaterialPropertyLib
{
/**
 * \class Phase
 * \brief This class defines material phases.
 * \details The Phase class consists of a vector of
 * components and an array of properties.
 */
class Phase
{
private:
    std::vector<std::unique_ptr<Component>> _components;
    PropertyArray _properties;
public:
    /// The Phase constructor is called with the optional
    /// phase name.
    Phase(boost::optional<std::string> const&);
    Phase();

    ~Phase(){};

    /// The method creating phase components based on config subtree.
    void createComponents(BaseLib::ConfigTree const&);
    /// The method creating phase properties based on config subtree.
    void createProperties(BaseLib::ConfigTree const&);
    /// The method initializing the phase properties.
    void createDefaultProperties(void);

    /// A simple get-function for a component. The argument refers to the
    /// Index of the component in the components vector.
    Component* component (std::size_t const&);
    /// A get-function for a property. The argument refers to the
    /// name of the property.
    Property* property(PropertyEnum const&);
    /// A get-function for retrieving the number of components in this phase
    std::size_t numberOfComponents(void);
    void resetPropertyUpdateStatus(void);
};

} //MaterialPropertyLib



#endif /* MATERIALLIB_MPL_MPPHASE_H_ */
