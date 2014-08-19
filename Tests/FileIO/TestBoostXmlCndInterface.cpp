/**
 * \file   TestBoostXmlCndInterface.cpp
 * \date   2014-02-07
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <boost/filesystem.hpp>

#include "gtest/gtest.h"

#include "BaseLib/BuildInfo.h"

// FileIO
#include "XmlIO/Qt/XmlCndInterface.h"
#include "XmlIO/Boost/BoostXmlCndInterface.h"

// OgsLib
#include "OGS/ProjectData.h"
#include "OGS/BoundaryCondition.h"
#include "OGS/FEMEnums.h"

TEST(FileIO, TestBoostXmlCndInterfaceUsingBoundaryCondition)
{
	// setup test data
	std::string geometry_name("GeometryForBC");
	const std::string bc_pnt_name("bc_pnt");
	GeoLib::Point *bc_pnt(new GeoLib::Point(0.0, 0.1, 0.2));

	std::vector<GeoLib::Point*> *pnts = new std::vector<GeoLib::Point*>;
	pnts->push_back(bc_pnt);

	std::map<std::string, std::size_t>* pnt_names
		= new std::map<std::string, std::size_t>;
	pnt_names->insert(std::pair<std::string, std::size_t>(bc_pnt_name, 0));

	ProjectData project_data;
	GeoLib::GEOObjects &geo_objects(*project_data.getGEOObjects());

	geo_objects.addPointVec(pnts, geometry_name, pnt_names);

	// fill BoundaryCondition data structure
	BoundaryCondition *bc(new BoundaryCondition(geometry_name));
	bc->initGeometricAttributes(geometry_name, GeoLib::GEOTYPE::POINT,
		bc_pnt_name, geo_objects);
	// set process info
	bc->setProcessType(FiniteElement::ProcessType::GROUNDWATER_FLOW);
	bc->setProcessPrimaryVariable(FiniteElement::PrimaryVariable::HEAD);
	// set distribution info
	bc->setProcessDistributionType(FiniteElement::DistributionType::CONSTANT);
	bc->setConstantDisValue(10.0);

	project_data.addCondition(bc);

	// write bc to file
	const std::string local_path("XmlCndInterfaceTestFile.cnd");
	std::string fname(BaseLib::BuildInfo::tests_tmp_path + local_path);

	FileIO::XmlCndInterface xml(project_data);
	xml.setNameForExport(geometry_name);
	int result_out = xml.writeToFile(fname);
	ASSERT_EQ(result_out, 1);

	// read bc from file using BoostXmlCndInterface
	FileIO::BoostXmlCndInterface cnd_interface(project_data);
	int result = cnd_interface.readFile(fname);

	//boost::filesystem::remove(fname);

	ASSERT_EQ(result, 1);
	// check the number of conditions in the vector
	std::vector<FEMCondition*> conds(project_data.getConditions());
	ASSERT_EQ(2u, conds.size());
	// compare the associated geometry and the names
	ASSERT_EQ(conds[0]->getAssociatedGeometryName().compare(
			conds[1]->getAssociatedGeometryName()
		), 0);
	ASSERT_EQ(conds[0]->getGeoName().compare(conds[1]->getGeoName()), 0);

	// fetch and compare geometries
	GeoLib::GeoObject * geo_obj_0 (
		const_cast<GeoLib::GeoObject*>(
			(geo_objects.getGeoObject(
				conds[0]->getAssociatedGeometryName(),
				GeoLib::GEOTYPE::POINT,
				conds[0]->getGeoName())
			)
		)
	);
	GeoLib::GeoObject * geo_obj_1 (
		const_cast<GeoLib::GeoObject*>(
			(geo_objects.getGeoObject(
				conds[1]->getAssociatedGeometryName(),
				GeoLib::GEOTYPE::POINT,
				conds[1]->getGeoName())
			)
		)
	);
	GeoLib::Point &pnt_0(* static_cast<GeoLib::Point *>(geo_obj_0));
	GeoLib::Point &pnt_1(* static_cast<GeoLib::Point *>(geo_obj_1));
	ASSERT_NEAR(pnt_0[0], pnt_1[0], std::numeric_limits<double>::epsilon());
	ASSERT_NEAR(pnt_0[1], pnt_1[1], std::numeric_limits<double>::epsilon());
	ASSERT_NEAR(pnt_0[2], pnt_1[2], std::numeric_limits<double>::epsilon());

	// compare process related information
	ASSERT_EQ(conds[0]->getProcessType(), conds[1]->getProcessType());
	ASSERT_EQ(conds[0]->getProcessPrimaryVariable(), conds[1]->getProcessPrimaryVariable());

	// compare distribution information
	ASSERT_EQ(conds[0]->getProcessDistributionType(),
		conds[1]->getProcessDistributionType());
	ASSERT_EQ(conds[0]->getDisValues().size(), conds[1]->getDisValues().size());
	ASSERT_NEAR(conds[0]->getDisValues()[0], conds[1]->getDisValues()[0], std::numeric_limits<double>::epsilon());
}
