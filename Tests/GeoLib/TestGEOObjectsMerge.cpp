/**
 * @file TestGEOObjectsMerge.cpp
 * @author Thomas Fischer
 * @date May 21, 2013
 * @brief 
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// google test include
#include "gtest/gtest.h"

// STL
#include <vector>
#include <map>
#include <string>

// BaseLib
#include "StringTools.h"

// GeoLib
#include "GEOObjects.h"

void createSetOfTestPointsAndAssociatedNames(GeoLib::GEOObjects & geo_objs, std::string &name, GeoLib::Point const& shift)
{
	std::vector<GeoLib::Point*> *pnts(new std::vector<GeoLib::Point*>);
	std::map<std::string, std::size_t>* pnt_name_map(new std::map< std::string, std::size_t>);

	const std::size_t pnts_per_edge(8);
	for (std::size_t k(0); k < pnts_per_edge; k++) {
		const std::size_t k_offset(k * pnts_per_edge * pnts_per_edge);
		for (std::size_t j(0); j < pnts_per_edge; j++) {
			const std::size_t offset(j * pnts_per_edge + k_offset);
			for (std::size_t i(0); i < pnts_per_edge; i++) {
				pnts->push_back(new GeoLib::Point(i+shift[0], j+shift[1], k+shift[2]));
				std::string pnt_name(
						name + "-" + BaseLib::number2str(i) + "-" + BaseLib::number2str(j) + "-"
								+ BaseLib::number2str(k));
				pnt_name_map->insert(std::pair< std::string, std::size_t>(pnt_name, i + offset));
			}
		}
	}

	geo_objs.addPointVec(pnts, name, pnt_name_map);
}

TEST(GeoLib, GEOObjectsMergePoints)
{
	GeoLib::GEOObjects geo_objs;
	std::vector<std::string> names;

	// *** insert set of points number 0
	GeoLib::Point shift (0.0,0.0,0.0);
	names.push_back("PointSet0");
	createSetOfTestPointsAndAssociatedNames(geo_objs, names[0], shift);

	// *** insert set of points number 1
	names.push_back("PointSet1");
	createSetOfTestPointsAndAssociatedNames(geo_objs, names[1], shift);

	// *** merge geometries
	std::string merged_geometries_name("MergedEqualPointSet");
	geo_objs.mergeGeometries(names, merged_geometries_name);

	GeoLib::PointVec const* merged_point_vec (geo_objs.getPointVecObj(merged_geometries_name));

	ASSERT_EQ(merged_point_vec->size(), 512);
	std::string test_name;
	merged_point_vec->getNameOfElementByID(0, test_name);
	ASSERT_EQ(test_name, "PointSet0-0-0-0");
	merged_point_vec->getNameOfElementByID(511, test_name);
	ASSERT_EQ(test_name, "PointSet0-7-7-7");

	// *** insert "shifted" set of points
	shift[0] += 8.0 * std::numeric_limits<double>::epsilon();
	names.push_back("ShiftedPointSet");
	createSetOfTestPointsAndAssociatedNames(geo_objs, names[2], shift);

	// *** merge PointSet0, PointSet1 and ShiftedPointSet
	merged_geometries_name = "MergedShiftedPointSet";
	geo_objs.mergeGeometries(names, merged_geometries_name);
	merged_point_vec = geo_objs.getPointVecObj(merged_geometries_name);

	ASSERT_EQ(merged_point_vec->size(), 1024);
	merged_point_vec->getNameOfElementByID(0, test_name);
	ASSERT_EQ(test_name, "PointSet0-0-0-0");
	merged_point_vec->getNameOfElementByID(511, test_name);
	ASSERT_EQ(test_name, "PointSet0-7-7-7");
	merged_point_vec->getNameOfElementByID(512, test_name);
	ASSERT_EQ(test_name, "ShiftedPointSet-0-0-0");
	merged_point_vec->getNameOfElementByID(1023, test_name);
	ASSERT_EQ(test_name, "ShiftedPointSet-7-7-7");

	std::size_t id;
	ASSERT_TRUE(merged_point_vec->getElementIDByName (test_name, id));
	ASSERT_EQ(id,1023);

	test_name = "PointSet1-0-0-0";
	ASSERT_FALSE(merged_point_vec->getElementIDByName (test_name, id));
}
