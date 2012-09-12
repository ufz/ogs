/*
 * PetrelInterface.h
 *
 *  Created on: Feb 16, 2010
 *      Author: fischeth
 */

#ifndef PETRELIO_H_
#define PETRELIO_H_

#include "GEOObjects.h"
#include <iostream>
#include <list>
#include <string>
#include <vector>

namespace FileIO
{
class PetrelInterface
{
public:
	PetrelInterface(std::list<std::string> &sfc_fnames,
	                std::list<std::string> &well_path_fnames,
	                std::string &unique_model_name,
	                GeoLib::GEOObjects* obj);
	virtual ~PetrelInterface();

private:
	void readPetrelSurface (std::istream &in);
	void readPetrelWellTrace (std::istream &in);
	void readPetrelWellTraceData (std::istream &in);
	std::string _unique_name;
	std::vector<GeoLib::Point*>* pnt_vec;
	std::vector<GeoLib::Point*>* well_vec;
	std::vector<GeoLib::Polyline*>* ply_vec;
	static const std::size_t MAX_COLS_PER_ROW = 256;
};
} // end namespace FileIO

#endif /* PETRELIO_H_ */
