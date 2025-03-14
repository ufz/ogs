/**
 * \file
 * \author Thomas Fischer
 * \date   2010-02-16
 * \brief  Definition of the PetrelInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

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
    PetrelInterface(std::list<std::string> const&sfc_fnames,
                    std::list<std::string> const&well_path_fnames,
                    std::string &unique_model_name,
                    GeoLib::GEOObjects* geo_obj);

    PetrelInterface(PetrelInterface const& other) = delete;
    PetrelInterface(PetrelInterface&& other) = delete;
    PetrelInterface& operator=(PetrelInterface const&) = delete;
    PetrelInterface& operator=(PetrelInterface&&) = delete;

private:
    void readPetrelSurfacePoints(std::istream& in);
    void readPetrelWellTrace (std::istream &in);
    void readPetrelWellTraceData (std::istream &in);
    std::string _unique_name;
    std::vector<GeoLib::Point*> pnt_vec;
    std::vector<GeoLib::Point*> well_vec;
};
} // end namespace FileIO
