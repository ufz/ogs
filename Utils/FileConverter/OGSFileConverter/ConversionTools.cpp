/**
 * \file ConversionTools.cpp
 * 2012/04/11 KR Initial implementation
 */

#include "ConversionTools.h"
#include "ProjectData.h"

#include <iostream>
#include <fstream>

// FEM Conditions
#include "BoundaryCondition.h"
#include "InitialCondition.h"
#include "SourceTerm.h"
#include "rf_bc_new.h"
#include "rf_ic_new.h"
#include "rf_st_new.h"
#include "FEMIO/BoundaryConditionIO.h"

// Qt
#include <QFileInfo>

void ConversionTools::getFEMConditionsFromASCIIFile(const QString &file_name, std::vector<FEMCondition*> &conditions)
{
	std::ifstream in(file_name.toStdString().data(), std::ios::in);
	if (!in.good())
	{
		std::cout << "Error in readASCIIConditionFile() - Could not open file." << std::endl;
		return;
	}

	QFileInfo fi(file_name);
	GeoLib::GEOObjects geo_objects;
	std::string geo_name(fi.baseName().toStdString() + ".gli");
	std::string file_path(fi.absolutePath().toStdString() + "/");
	FEMCondition::CondType type(FEMCondition::UNSPECIFIED);
	std::string cond_tag("");
	if (fi.suffix().toLower() == "bc")
	{
		cond_tag = "#BOUNDARY_CONDITION";
		type = FEMCondition::BOUNDARY_CONDITION;
	}
	else if (fi.suffix().toLower() == "ic")
	{
		cond_tag = "#INITIAL_CONDITION";
		type = FEMCondition::INITIAL_CONDITION;
	}
	else if (fi.suffix().toLower() == "st")
	{
		cond_tag = "#SOURCE_TERM";
		type = FEMCondition::SOURCE_TERM;
	}

	std::cout << "Reading " << fi.fileName().toStdString() << "..." << std::endl;
	while (!in.eof())
	{
		char buffer[256];
		in.getline(buffer, 256);
		std::string line(buffer);
		if (line.find("#STOP") != std::string::npos)
			return;
		if (line.find(cond_tag) != std::string::npos)
		{
			std::ios::pos_type position = in.tellg();
			if (type == FEMCondition::BOUNDARY_CONDITION)     position = readBoundaryCondition(conditions, in, file_path, geo_objects, geo_name);
			else if (type == FEMCondition::INITIAL_CONDITION) position = readInitialCondition(conditions, in, file_path, geo_objects, geo_name);
			else if (type == FEMCondition::SOURCE_TERM)       position = readSourceTerm(conditions, in, file_path, geo_objects, geo_name);
			in.seekg(position, std::ios::beg);
		}
	}
}

std::ios::pos_type ConversionTools::readBoundaryCondition(std::vector<FEMCondition*> &conditions, std::ifstream &in, const std::string &file_path, const GeoLib::GEOObjects &geo_objects, const std::string &geo_name)
{
	Q_UNUSED(file_path);
	CBoundaryCondition* bc(new CBoundaryCondition());
	bool valid;
	std::ios::pos_type position = bc->Read(&in, geo_objects, geo_name, valid);
	conditions.push_back(new BoundaryCondition(*bc, geo_name));
	delete bc;
	return position;
}

std::ios::pos_type ConversionTools::readInitialCondition(std::vector<FEMCondition*> &conditions, std::ifstream &in, const std::string &file_path, const GeoLib::GEOObjects &geo_objects, const std::string &geo_name)
{
	Q_UNUSED(file_path);
	CInitialCondition* ic = new CInitialCondition();
	std::ios::pos_type position = ic->Read(&in, geo_objects, geo_name);
	conditions.push_back(new InitialCondition(*ic, geo_name));
	delete ic;
	return position;
}

std::ios::pos_type ConversionTools::readSourceTerm(std::vector<FEMCondition*> &conditions, std::ifstream &in, const std::string &file_path, const GeoLib::GEOObjects &geo_objects, const std::string &geo_name)
{
	CSourceTerm* st(new CSourceTerm());
	std::ios::pos_type position = st->Read(&in, geo_objects, geo_name);
	conditions.push_back(new SourceTerm(*st, geo_name));

	if (st->getProcessDistributionType() == FiniteElement::DIRECT)
		conditions[conditions.size()-1]->setDisValues( ConversionTools::getDirectNodeValues(file_path + st->fname) );

	delete st;
	return position;
}

std::vector< std::pair<size_t, double> > ConversionTools::getDirectNodeValues(std::string file_name)
{
	std::vector< std::pair<size_t, double> > node_values;
	SourceTerm::getDirectNodeValues(file_name, node_values);
	return node_values;
}

int ConversionTools::writeDirectValues(const FEMCondition &condition, const std::string &direct_value_file)
{
	std::ofstream out(direct_value_file.c_str());
	if (!out.good())
	{
		std::cout << "Error in writeDirectValues() - Could not open file." << std::endl;
		return 0;
	}

	const std::vector<size_t> dis_nodes = condition.getDisNodes();
	const std::vector<double> dis_values = condition.getDisValues();
	const size_t nValues(dis_nodes.size());

	if (nValues != dis_values.size() || nValues==0)
		return 0;

	for (size_t i=0; i<nValues; i++)
		out << dis_nodes[i] << "\t" << dis_values[i] << std::endl;

	return 1;
}

