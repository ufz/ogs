/**
 * \file ProjectData.cpp
 * 25/08/2010 KR Initial implementation
 */

#include "ProjectData.h"
#include "StringTools.h"


ProjectData::ProjectData()
//: _geoObjects ()
{}

ProjectData::~ProjectData()
{
	delete _geoObjects;
	for (std::map<std::string, MeshLib::CFEMesh*>::iterator it = _msh_vec.begin(); it != _msh_vec.end(); ++it)
	{
		delete it->second;
	}
	size_t nCond (_cond_vec.size());
	for (size_t i=0; i<nCond; i++)
	{
		delete _cond_vec[i];
	}
}

void ProjectData::addMesh(MeshLib::CFEMesh* mesh, std::string &name)
{
	isUniqueMeshName(name);
	_msh_vec[name] = mesh;
};

const MeshLib::CFEMesh* ProjectData::getMesh(const std::string &name) const
{
	return _msh_vec.find(name)->second;
}

bool ProjectData::removeMesh(const std::string &name)
{
	delete _msh_vec[name];
	size_t result = _msh_vec.erase(name);
	return (result>0);
}

void ProjectData::addCondition(FEMCondition* cond)
{
	_cond_vec.push_back(cond);
};

void ProjectData::addConditions(std::vector<FEMCondition*> conds)
{
	for (size_t i=0; i<conds.size(); i++)
		_cond_vec.push_back(conds[i]);
};

const FEMCondition* ProjectData::getCondition(const std::string &geo_name, GEOLIB::GEOTYPE type, const std::string &cond_name) const
{
	for (std::vector<FEMCondition*>::const_iterator it = _cond_vec.begin(); it != _cond_vec.end(); ++it)
	{
		if ((*it)->getAssociatedGeometryName().compare(geo_name) == 0) 
		{
			if ( ((*it)->getGeoName().compare(cond_name)==0) && ((*it)->getGeoType()==type) )
				return *it;
		}
	}
	std::cout << "Error in ProjectData::getCondition() - No condition found with name \"" << cond_name << "\"..." << std::endl;
	return NULL;
}

const std::vector<FEMCondition*> ProjectData::getConditions(const std::string &geo_name, FEMCondition::CondType type) const
{
	std::vector<FEMCondition*> conds;
	for (std::vector<FEMCondition*>::const_iterator it = _cond_vec.begin(); it != _cond_vec.end(); ++it)
	{
		if ((*it)->getAssociatedGeometryName().compare(geo_name) == 0)
		{
			if ( (type == FEMCondition::UNSPECIFIED) || ((*it)->getCondType() == type) )
				conds.push_back(*it);
		}
	}
	return conds;
}

bool ProjectData::removeCondition(const std::string &geo_name, GEOLIB::GEOTYPE type, const std::string &cond_name)
{
	for (std::vector<FEMCondition*>::iterator it = _cond_vec.begin(); it != _cond_vec.end(); ++it)
	{
		if ((*it)->getAssociatedGeometryName().compare(geo_name) == 0) 
		{
			if ( ((*it)->getGeoName().compare(cond_name)==0) && ((*it)->getGeoType()==type) )
			{
				delete *it;
				_cond_vec.erase(it);
				return true;
			}
		}
	}
	std::cout << "Error in ProjectData::getCondition() - No condition found with name \"" << cond_name << "\"..." << std::endl;
	return false;
}

void ProjectData::removeConditions(const std::string &geo_name, FEMCondition::CondType type)
{
	for (std::vector<FEMCondition*>::iterator it = _cond_vec.begin(); it != _cond_vec.end();)
	{
		if ( ((*it)->getAssociatedGeometryName().compare(geo_name) == 0) 
			&& ( (type == FEMCondition::UNSPECIFIED) || ((*it)->getCondType() == type) ))
		{
			delete *it;
			it = _cond_vec.erase(it);
		}
		else ++it;
	}
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
		if (count>1) cpName = cpName + "-" + number2str(count);

		for (std::map<std::string, MeshLib::CFEMesh*>::iterator it = _msh_vec.begin(); it != _msh_vec.end(); ++it)
		{
			if ( cpName.compare(it->first) == 0 ) isUnique = false;
		}
	}

	// At this point cpName is a unique name and isUnique is true.
	// If cpName is not the original name, "name" is changed and isUnique is set to false,
	// indicating that a vector with the original name already exists.
	if (count>1)
	{
		isUnique = false;
		name = cpName;
	}
	return isUnique;
}
