/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FEMCondition.h
 *
 * Created on 2010-11-25 by Karsten Rink
 *
 */

#ifndef FEMCONDITION_H
#define FEMCONDITION_H

#include "DistributionInfo.h"
#include "GeoInfo.h"
#include "GeoObject.h"
#include "Point.h"
#include "ProcessInfo.h"

#include <vector>

class CBoundaryCondition;
class CInitialCondition;
class CSourceTerm;

//class GEOObjects;
class GridAdapter;
#include "GEOObjects.h"

/**
 * \brief Adapter class for handling FEM Conditions in the user Interface
 */
class FEMCondition : public ProcessInfo, public GeoInfo, public DistributionInfo
{
public:
	/// Specifier for types of FEM Conditions
	enum CondType {
		UNSPECIFIED        = 0,
		BOUNDARY_CONDITION = 1,
		INITIAL_CONDITION  = 2,
		SOURCE_TERM        = 3
	};

	FEMCondition(const std::string &geometry_name, CondType t = UNSPECIFIED);
	FEMCondition(const std::string &geometry_name,
					FiniteElement::ProcessType pt = FiniteElement::INVALID_PROCESS,
					FiniteElement::PrimaryVariable pv =	FiniteElement::INVALID_PV,
					GeoLib::GEOTYPE gt = GeoLib::INVALID,
					const std::string &gn = "[unspecified]",
					FiniteElement::DistributionType dt = FiniteElement::INVALID_DIS_TYPE, CondType ct = UNSPECIFIED);
	FEMCondition(const FEMCondition &cond, CondType t);

	~FEMCondition() {}

	/// Returns the type of the FEM Condition (i.e. BC, IC or ST)
	CondType getCondType() const { return _type; }

	/// Returns the value(s) for the distribution
	const std::vector<std::size_t> getDisNodes() const { return _disNodes; }

	/// Returns the value(s) for the distribution
	const std::vector<double> getDisValues() const { return _disValues; }

	/// Returns the name of the geo-object the condition is assigned to. This object is part of the associated geometry.
	const std::string& getGeoName() const { return _geoName; }

	/// Returns the name of the associated geometry.
	const std::string& getAssociatedGeometryName() const { return _associated_geometry; }

	/// Convenience method for setting a single value specifying the distribution.
	void setConstantDisValue(double disValue) {_disValues.clear(); _disValues.push_back(disValue); }

	/// Sets a vector of values specifying the distribution.
	void setDisValues(const std::vector<std::size_t> &disNodes, const std::vector<double> &disValues)
	{
		_disNodes = disNodes;
		_disValues = disValues;
	}

	/// Sets a vector of values specifying the distribution.
	/// The first value specifies the point id, the second the value for that point.
	void setDisValues(const std::vector< std::pair<std::size_t, double> > &dis_values);

	/// Removes all distribution values
	void clearDisValues() { _disValues.resize(0); };

	/// Sets the name of the geo-object the condition is assigned to.
	void setGeoName(std::string geoName) { _geoName = geoName; }

	/// Returns the type of the FEM condition as a string.
	static std::string condTypeToString(CondType type);

protected:
	CondType _type;
	std::string _geoName;
	std::vector<std::size_t> _disNodes;
	std::vector<double> _disValues;
	std::string _associated_geometry;
};

#endif //FEMCONDITION_H
