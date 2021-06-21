/**
 * \file
 * \author Karsten Rink
 * \date   2013-03-18
 * \brief  Definition of the StationBorehole class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Station.h"

#include <list>
#include <string>
#include <vector>

namespace GeoLib
{

/**
 * \brief A borehole as a geometric object.
 *
 * A borehole inherits Station but has some additional information such as a date, a borehole profile, etc.
 */
class StationBorehole : public Station
{
public:
    /** constructor initialises the borehole with the given coordinates */
    explicit StationBorehole(double x = 0.0,
                             double y = 0.0,
                             double z = 0.0,
                             double const depth = 0.0,
                             const std::string& name = "");
    ~StationBorehole() override;

    /// Creates a StationBorehole-object from a string (assuming the string has the right format)
    static StationBorehole* createStation(const std::string &line);

    /// Creates a new borehole object based on the given parameters.
    static StationBorehole* createStation(const std::string &name,
                                          double x,
                                          double y,
                                          double z,
                                          double depth,
                                          const std::string &date = "");

    // Returns the depth of the borehole
    double getDepth() const { return _depth; }

    /// Returns the date entry for the borehole
    double getDate() const { return _date; }

    /// Returns a reference to a vector of Points representing the stratigraphy of the borehole (incl. the station-point itself)
    const std::vector<Point*> &getProfile() const { return _profilePntVec; }

    /// Returns a reference to a vector of soil names for the stratigraphy of the borehole
    const std::vector<std::string> &getSoilNames() const { return _soilName; }

    /// Sets the depth of the borehole
    void setDepth( double depth ) { _depth = depth; }

    /// Add a soil layer to the boreholes stratigraphy.
    void addSoilLayer ( double thickness, const std::string &soil_name);

    /**
     * Add a soil layer to the boreholes stratigraphy.
     * Note: The given coordinates always mark THE END of the soil layer. The reason behind this is
     * that the beginning of the first layer is identical with the position of the borehole. For each
     * layer following the beginning is already given by the end of the last layer. This also saves
     * a separate entry in the profile vector for the end of the borehole which in the given notation
     * is just the coordinate given for the last soil layer (i.e. the end of that layer).
     */
    void addSoilLayer ( double x, double y, double z, const std::string &soil_name);

private:
    //long profile_type;
    //std::vector<long> _soilType;
    double _depth{0};  // depth of the borehole
    int _date{0};      // date when the borehole has been drilled

    /// Contains the names for all the soil layers
    std::vector<std::string> _soilName;

    /// Contains the points for the lower boundaries of all layers
    std::vector<Point*> _profilePntVec;
};

bool isBorehole(GeoLib::Point const* pnt);

}  // namespace GeoLib
