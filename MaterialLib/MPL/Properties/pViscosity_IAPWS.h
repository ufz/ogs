/**
 * \author Norbert Grunwald
 * \date   Sep 21, 2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#ifndef MATERIALLIB_MPL_PROPERTIES_PVISCOSITY_IAPWS_H_
#define MATERIALLIB_MPL_PROPERTIES_PVISCOSITY_IAPWS_H_

#include "../mpProperty.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;

/**
 *
 *  \brief J.R. Cooper et. al: The International Association for the Properties
 *        of Water and Steam. Berlin, Germany, September 2008.
 */
class ViscosityWaterIAPWS final : public Property
{
private:
    Component* _component;
    std::array<double, 4> _H;
    std::array<std::array<double, 7>, 6> _h;

public:
    ViscosityWaterIAPWS(Medium*);
    ViscosityWaterIAPWS(Phase*);
    ViscosityWaterIAPWS(Component*);

    /// This method overrides the base class implementation and
    /// actually computes and sets the property _value.
    PropertyDataType value(VariableArray const&) override;
};

}  // MaterialPropertyLib

#endif /* MATERIALLIB_MPL_PROPERTIES_PVISCOSITY_IAPWS_H_ */
