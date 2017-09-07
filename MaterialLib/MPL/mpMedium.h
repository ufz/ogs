/**
 * \author Norbert Grunwald
 * \date   07.09.2017
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "BaseLib/ConfigTree.h"
#include "mpPhase.h"
#include "mpProperty.h"

#include <memory>
#include <vector>

namespace MaterialPropertyLib
{
/// This class is for material objects on the Medium scale.
///
/// A Medium consists of an arbitrarily long vector of phases and an array of
/// properties.
class Medium
{
private:
    /// The vector that holds the phases.
    std::vector<std::unique_ptr<Phase>> _phases;
    /// The array that holds the medium properties.
    PropertyArray _properties;

public:
    /// This constructor parses the "phases" and "properties" subtrees of the
    /// config tree and calls create methods for the phase vector and the
    /// properties array. Medium properties are optional. If not defined,
    /// default properties are assigned.
    Medium(BaseLib::ConfigTree const& config);

    /// A method that parses the phase details and stores them in the
    /// private _phases member.
    ///
    /// This method creates the phases of the medium. Unlike a medium, a phase
    /// may have a name. However, this is silly at the moment since this name
    /// still has no effect (except of some benefits in regard of readability).
    /// Phase components are required (a phase consists of at least one
    /// component). Phase properties are optional. If not given, default
    /// properties are assigned. These default properties average the component
    /// properties, weighted by mole fraction.
    void createPhases(BaseLib::ConfigTree const& config);

    /// A method that parses the medium property details and stores
    /// them in the private _properties member.
    ///
    /// This method creates the properties of the Medium as defined in the
    /// prj-file. Only specified properties overwrite the default properties.
    void createProperties(BaseLib::ConfigTree const& config);

    /// A method that creates the default properties of the medium.  Currently,
    /// these defaults is the volume fraction average.
    ///
    /// Most properties are fine with the volume fraction average, but
    /// special-defaults are allowed as well...
    void createDefaultProperties();

    /// A get-function for a particular phase. The ul argument specifies the
    /// index within the _phases vector.
    Phase& phase(std::size_t index) const;
    /// A get-function for a particular phase by phase name.
    Phase& phase(std::string const& phase_name) const;
    /// A get-function for a property. The argument refers to the name of the
    /// property.
    Property& property(PropertyEnum const& p) const;

    /// A simple get-function for retrieving the number of phases the medium
    /// consists of.
    std::size_t numberOfPhases() const;

};
}  // namespace MaterialPropertyLib
