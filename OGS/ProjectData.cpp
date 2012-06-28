/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
 * \file ProjectData.cpp
 *
 * Created on 2010-08-25 by Karsten Rink
 */

#include "ProjectData.h"
#include "StringTools.h"

ProjectData::ProjectData()
: _geoObjects (NULL)
{}

ProjectData::~ProjectData()
{
	delete _geoObjects;
	for (std::map<std::string, MeshLib::Mesh*>::iterator it = _msh_vec.begin();
	     it != _msh_vec.end(); ++it)
		delete it->second;
	size_t nCond (_cond_vec.size());
	for (size_t i = 0; i < nCond; i++)
		delete _cond_vec[i];
}

void ProjectData::addMesh(MeshLib::Mesh* mesh, std::string &name)
{
	isUniqueMeshName(name);
	_msh_vec[name] = mesh;
}

const MeshLib::Mesh* ProjectData::getMesh(const std::string &name) const
{
	return _msh_vec.find(name)->second;
}

bool ProjectData::removeMesh(const std::string &name)
{
	delete _msh_vec[name];
	size_t result = _msh_vec.erase(name);
	return result > 0;
}

bool ProjectData::meshExists(const std::string &name)
{
	if (_msh_vec.count(name)>0) return true;
	return false;
}

void ProjectData::addProcess(ProcessInfo* pcs)
{
	for (std::vector<ProcessInfo*>::const_iterator it = _pcs_vec.begin(); it != _pcs_vec.end(); ++it)
		if ((*it)->getProcessType() == pcs->getProcessType())
		{
			std::cout << "Error in ProjectData::addProcess() - "
			  << FiniteElement::convertProcessTypeToString(pcs->getProcessType()) << " process already exists." << std::endl;
		}

	_pcs_vec.push_back(pcs);
}

const ProcessInfo* ProjectData::getProcess(FiniteElement::ProcessType type) const
{
	for (std::vector<ProcessInfo*>::const_iterator it = _pcs_vec.begin(); it != _pcs_vec.end(); ++it)
		if ((*it)->getProcessType() == type)
			return *it;

	std::cout << "Error in ProjectData::getProcess() - No "
			  << FiniteElement::convertProcessTypeToString(type) << " process found..." << std::endl;
	return NULL;
}

bool ProjectData::removeProcess(FiniteElement::ProcessType type)
{
	for (std::vector<ProcessInfo*>::iterator it = _pcs_vec.begin(); it != _pcs_vec.end(); ++it)
		if ((*it)->getProcessType() == type)
		{

			delete *it;
			_pcs_vec.erase(it);
			return true;
		}

	std::cout << "Error in ProjectData::removeProcess() - No "
			  << FiniteElement::convertProcessTypeToString(type) << " process found..." << std::endl;
	return false;
}

void ProjectData::addCondition(FEMCondition* cond)
{
	_cond_vec.push_back(cond);
}

void ProjectData::addConditions(std::vector<FEMCondition*> conds)
{
	for (size_t i = 0; i < conds.size(); i++)
		_cond_vec.push_back(conds[i]);
}

const FEMCondition* ProjectData::getCondition(const std::string &geo_name,
                                              GeoLib::GEOTYPE type,
                                              const std::string &cond_name) const
{
	for (std::vector<FEMCondition*>::const_iterator it = _cond_vec.begin(); it != _cond_vec.end(); ++it)
		if ((*it)->getAssociatedGeometryName().compare(geo_name) == 0)
			if ( ((*it)->getGeoName().compare(cond_name) == 0) &&
			     ((*it)->getGeoType() == type) )
				return *it;

	std::cout << "Error in ProjectData::getCondition() - No condition found with name \"" <<
	cond_name << "\"..." << std::endl;
	return NULL;
}

const std::vector<FEMCondition*> ProjectData::getConditions(FiniteElement::ProcessType pcs_type,
															std::string geo_name,
                                                            FEMCondition::CondType cond_type) const
{
	// if all
	if (pcs_type == FiniteElement::INVALID_PROCESS && geo_name.empty() && cond_type == FEMCondition::UNSPECIFIED)
		return _cond_vec;

	// else: filter according to parameters
	std::vector<FEMCondition*> conds;
	for (std::vector<FEMCondition*>::const_iterator it = _cond_vec.begin(); it != _cond_vec.end(); ++it)
	{
		if ( ((pcs_type == FiniteElement::INVALID_PROCESS) || (pcs_type == ((*it)->getProcessType()))) &&
			 ((geo_name.empty() || ((*it)->getAssociatedGeometryName().compare(geo_name) == 0))) &&
			 ((cond_type == FEMCondition::UNSPECIFIED) || ((*it)->getCondType() == cond_type)) )
				conds.push_back(*it);
	}

	return conds;
}

void ProjectData::removeConditions(FiniteElement::ProcessType pcs_type, std::string geo_name, FEMCondition::CondType cond_type)
{
	// if all
	if (pcs_type == FiniteElement::INVALID_PROCESS && geo_name.empty() && cond_type == FEMCondition::UNSPECIFIED)
	{
		for (size_t i=0; i<_cond_vec.size(); i++) delete _cond_vec[i];
		return;
	}

	// else: filter according to parameters
	for (std::vector<FEMCondition*>::iterator it = _cond_vec.begin(); it != _cond_vec.end(); )
	{
		if ( ((pcs_type == FiniteElement::INVALID_PROCESS) || (pcs_type == ((*it)->getProcessType()))) &&
			 ((geo_name.empty() || ((*it)->getAssociatedGeometryName().compare(geo_name) == 0))) &&
			 ((cond_type == FEMCondition::UNSPECIFIED) || ((*it)->getCondType() == cond_type)) )
		{
			delete *it;
			it = _cond_vec.erase(it);
		}
		else
			++it;
	}
}

bool ProjectData::removeCondition(const std::string &geo_name,
                                  GeoLib::GEOTYPE type,
                                  const std::string &cond_name)
{
	for (std::vector<FEMCondition*>::iterator it = _cond_vec.begin(); it != _cond_vec.end(); ++it)
	{
		if ((*it)->getAssociatedGeometryName().compare(geo_name) == 0)
			if ( ((*it)->getGeoName().compare(cond_name) == 0) &&
			     ((*it)->getGeoType() == type) )
			{
				delete *it;
				_cond_vec.erase(it);
				return true;
			}
	}
	std::cout << "Error in ProjectData::removeCondition() - No condition found with name \"" <<
	cond_name << "\"..." << std::endl;
	return false;
}

bool ProjectData::isUniqueMeshName(std::string &name)
{
	int count(0);
	bool isUnique(false);
	std::string cpName;

	while (!isUnique)
	{
		isUnique = true;
		cpName = name;

		count++;
		// If the original name already exists we start to add numbers to name for
		// as long as it takes to make the name unique.
		if (count > 1)
			cpName = cpName + "-" + number2str(count);

		for (std::map<std::string, MeshLib::Mesh*>::iterator it = _msh_vec.begin();
		     it != _msh_vec.end(); ++it)
			if ( cpName.compare(it->first) == 0 )
				isUnique = false;
	}

	// At this point cpName is a unique name and isUnique is true.
	// If cpName is not the original name, "name" is changed and isUnique is set to false,
	// indicating that a vector with the original name already exists.
	if (count > 1)
	{
		isUnique = false;
		name = cpName;
	}
	return isUnique;
}
