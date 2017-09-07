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

#include "mpProperty.h"
#include "Properties/properties.h"


namespace MaterialPropertyLib
{

Property::Property(){};

Property* newProperty(BaseLib::ConfigTree const& config)
{
        /// Parsing the property type:
        auto const property_type = config.getConfigParameter<std::string>("type");

        /// If (and only if) the given property type is 'constant', a
        /// corresponding value is needed.
        if (boost::iequals(property_type , "constant"))
        {
            auto const property_value = config.getConfigParameter<double>("value");
            return new Constant();
        }
        OGS_FATAL("The specified component property type \"%s\" was not "
                "recognized", property_type.c_str());
        return nullptr;
}


} // MaterialPropertyLib
