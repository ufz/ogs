/**
 * \author Norbert Grunwald
 * \date   18.09.2017
 * \brief  
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#ifndef MATERIALLIB_MPL_PROPERTIES_PPENGROBINSON_H_
#define MATERIALLIB_MPL_PROPERTIES_PPENGROBINSON_H_

#include "../mpProperty.h"

namespace MaterialPropertyLib
{

class Medium;
class Phase;
class Component;
/**
 * \class AverageVolumeFraction
 * \brief A function averaging a property by volume fraction
 * \details This property is usually a medium property, it
 * computes the average of individual phase properties
 * weighted by volume fraction.
 */
class PengRobinson : public Property
{
private:
	/// A pointer to the phase object.
	Phase* _phase;
	/// A pointer to the component object.
	Component* _component;

	std::array<double, 2> _a {{0}}; // cohesion pressure
	std::array<double, 2> _b {{0}}; // co-volume
	std::array<double, 2> _k {{0}}; // co-volume


public:
	/// Constructor passing a pointer to the medium.
	PengRobinson(Medium*);
	/// Constructor passing a pointer to a phase.
	PengRobinson(Phase*);
	/// Constructor passing a pointer to a component.
	PengRobinson(Component*);

	/// This method overrides the base class implementation and
	/// actually computes and sets the property _value.
	PropertyDataType value (VariableArray const&) override;
};

/// characteristic constant kappa as function of acentric factor omega
double kappa(const double);
/// dimensionless function alpha of reduced temperature
double alpha(const double, const double, const double);
double cohesionPressure (const double Tc, const double pc);
double coVolume (const double Tc, const double pc);

}  // MaterialPropertyLib




#endif /* MATERIALLIB_MPL_PROPERTIES_PPENGROBINSON_H_ */
