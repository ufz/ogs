/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TINInterface.h"

#include <fstream>
#include <limits>

#include <logog/include/logog.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "GeoLib/AnalyticalGeometry.h"
#include "GeoLib/Surface.h"
#include "GeoLib/Triangle.h"

#include "MathLib/GeometricBasics.h"

namespace GeoLib
{
namespace IO
{

GeoLib::Surface* TINInterface::readTIN(std::string const& fname,
    GeoLib::PointVec &pnt_vec,
    std::vector<std::string>* errors)
{
    // open file
    std::ifstream in(fname.c_str());
    if (!in) {
        WARN("readTIN(): could not open stream from %s.", fname.c_str());
        if (errors)
            errors->push_back ("readTINFile error opening stream from " + fname);
        return nullptr;
    }

    auto* sfc = new GeoLib::Surface(*(pnt_vec.getVector()));
    std::size_t id;
    MathLib::Point3d p0, p1, p2;
    std::string line;
    while (std::getline(in, line).good())
    {
        // allow empty lines
        if (line.empty())
            continue;

        // parse line
        std::stringstream input(line);
        // read id
        if (!(input >> id)) {
            in.close();
            delete sfc;
            return nullptr;
        }
        // read first point
        if (!(input >> p0[0] >> p0[1] >> p0[2])) {
            ERR("Could not read coords of 1st point of triangle %d.", id);
            if (errors)
                errors->push_back (std::string("readTIN error: ") +
                    std::string("Could not read coords of 1st point in triangle ") +
                    std::to_string(id));
            in.close();
            delete sfc;
            return nullptr;
        }
        // read second point
        if (!(input >> p1[0] >> p1[1] >> p1[2])) {
            ERR("Could not read coords of 2nd point of triangle %d.", id);
            if (errors)
                errors->push_back (std::string("readTIN error: ") +
                    std::string("Could not read coords of 2nd point in triangle ") +
                    std::to_string(id));
            in.close();
            delete sfc;
            return nullptr;
        }
        // read third point
        if (!(input >> p2[0] >> p2[1] >> p2[2])) {
            ERR("Could not read coords of 3rd point of triangle %d.", id);
            if (errors)
                errors->push_back (std::string("readTIN error: ") +
                    std::string("Could not read coords of 3rd point in triangle ") +
                    std::to_string(id));
            in.close();
            delete sfc;
            return nullptr;
        }

        // check area of triangle
        double const d_eps(std::numeric_limits<double>::epsilon());
        if (MathLib::calcTriangleArea(p0, p1, p2) < d_eps) {
            ERR("readTIN: Triangle %d has zero area.", id);
            if (errors)
                errors->push_back (std::string("readTIN: Triangle ")
                    + std::to_string(id) + std::string(" has zero area."));
            delete sfc;
            return nullptr;
        }

        // determine size pnt_vec to insert the correct ids
        std::size_t const s(pnt_vec.getVector()->size());

        std::size_t const pnt_pos_0(pnt_vec.push_back(new GeoLib::Point(p0, s)));
        std::size_t const pnt_pos_1(pnt_vec.push_back(new GeoLib::Point(p1, s+1)));
        std::size_t const pnt_pos_2(pnt_vec.push_back(new GeoLib::Point(p2, s+2)));
        // create new Triangle
        if (pnt_pos_0 != std::numeric_limits<std::size_t>::max() &&
            pnt_pos_1 != std::numeric_limits<std::size_t>::max() &&
            pnt_pos_1 != std::numeric_limits<std::size_t>::max()) {
            sfc->addTriangle(pnt_pos_0, pnt_pos_1, pnt_pos_2);
        }
    }

    if (sfc->getNumberOfTriangles() == 0) {
        WARN("readTIN(): No triangle found.", fname.c_str());
        if (errors)
            errors->push_back ("readTIN error because of no triangle found");
        delete sfc;
        return nullptr;
    }

    return sfc;
}

void TINInterface::writeSurfaceAsTIN(GeoLib::Surface const& surface, std::string const& file_name)
{
    std::ofstream os (file_name.c_str());
    if (!os) {
        WARN("writeSurfaceAsTIN(): could not open stream to %s.", file_name.c_str());
        return;
    }
    os.precision(std::numeric_limits<double>::digits10);
    const std::size_t n_tris (surface.getNumberOfTriangles());
    for (std::size_t l(0); l < n_tris; l++) {
        GeoLib::Triangle const& tri (*(surface[l]));
        os << l << " " << *(tri.getPoint(0)) << " " << *(tri.getPoint(1)) << " " << *(tri.getPoint(2)) << "\n";
    }
    os.close();
}
} // end namespace IO
} // end namespace GeoLib
