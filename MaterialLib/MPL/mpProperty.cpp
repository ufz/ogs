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

PropertyDataType Property::value()
{
    return _value;
}

Property* newProperty(BaseLib::ConfigTree const& config)
{
        /// Parsing the property type:
        auto const property_type = config.getConfigParameter<std::string>("type");

        /// If (and only if) the given property type is 'constant', a
        /// corresponding value is needed.
        if (boost::iequals(property_type , "constant"))
        {
            /// \todo: property_value is, at the moment, restricted to
            /// scalars, although the property class and its derivatives
            /// excepts the datatype PropertyDataType. If necessary, this
            /// method can be adapted to parse multiple value components,
            /// e.g. by
            /// <value> 0.0 0.0 -9.81</value>

            auto const property_value = config.getConfigParameter<double>("value");
            return new Constant(property_value);
        }
        OGS_FATAL("The specified component property type \"%s\" was not "
                "recognized", property_type.c_str());
        return nullptr; // to avoid 'no return' warnings
}


} // MaterialPropertyLib
