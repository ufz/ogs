/**
 * \file
 * \author Thomas Fischer
 * \date   2010-02-16
 * \brief  Definition of the PetrelInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * @file PetrelInterface.h
 * @date 2010-02-16
 * @author Thomas Fischer
 */

#ifndef PETRELIO_H_
#define PETRELIO_H_

#include <iosfwd>
#include <list>
#include <string>
#include <vector>

namespace GeoLib
{
    class GEOObjects;
    class Point;
    class Polyline;
}

namespace FileIO
{
class PetrelInterface final
{
public:
    PetrelInterface(std::list<std::string> &sfc_fnames,
                    std::list<std::string> &well_path_fnames,
                    std::string &unique_model_name,
                    GeoLib::GEOObjects* obj);

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
