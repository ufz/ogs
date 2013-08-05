/**
 * \file
 * \author Karsten Rink
 * \date   2013-03-18
 * \brief  Definition of the StationBorehole class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#ifndef GEO_STATIONBOREHOLE_H
#define GEO_STATIONBOREHOLE_H

#include "Station.h"

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
	StationBorehole(double x = 0.0, double y = 0.0, double z = 0.0, const std::string &name = "");
	~StationBorehole(void);

	/// Creates a StationBorehole-object from a string (assuming the string has the right format)
	static StationBorehole* createStation(const std::string &line);

	/// Creates a new borehole object based on the given parameters.
	static StationBorehole* createStation(const std::string &name,
	                                      double x,
	                                      double y,
	                                      double z,
	                                      double depth,
	                                      std::string date = "");

	/// Adds a stratigraphy to a borehole given a vector of points of length "n" and a vector of soil names of length "n-1".
	int addStratigraphy(const std::vector<Point*> &profile, const std::vector<std::string> &soil_names);

	/// Reads the stratigraphy for a specified station from a file
	static int addStratigraphy(const std::string &path, StationBorehole* borehole);

	/**
	 * \brief Reads all stratigraphy information from a file in one go.
	 *
	 * Reads all stratigraphy information from a file in one go.
	 * Be very careful when using this method -- it is pretty fast but it checks nothing and just
	 * assumes that everything is in right order and will work out fine!
	 */
	static int addStratigraphies(const std::string &path, std::vector<Point*>* boreholes);

	/// Finds the given string in the vector of soil-names
	int find(const std::string &str);

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
	 * a seperate entry in the profile vector for the end of the borehole which in the given notation
	 * is just the coordinate given for the last soil layer (i.e. the end of that layer).
	 */
	void addSoilLayer ( double x, double y, double z, const std::string &soil_name);

protected:
	/// Returns the depth of this borehole. Please see the documentation for Station::getX for details concerning the syntax.
	static double getDepth(void* stnObject)  { StationBorehole* stn =
		                                           (StationBorehole*)stnObject;
		                                   return stn->_depth; }
	/// Returns the date this borehole has been drilled. Please see the documentation for Station::getX for details concerning the syntax.
	static double getDate(void* stnObject)  { StationBorehole* stn = (StationBorehole*)stnObject; return stn->_date; }
	/// Sets the depth of this borehole. Please see the documentation for Station::getX for details concerning the syntax.
	static void setDepth(void* stnObject, double val) { StationBorehole* stn =
		                                                    (StationBorehole*)stnObject;
		                                            stn->_depth = val; }
	/// Sets the date when this borehole has been drilled. Please see the documentation for Station::getX for details concerning the syntax.
	static void setDate(void* stnObject, double val) { StationBorehole* stn = (StationBorehole*)stnObject; stn->_date = static_cast<int>(val); }

private:
	/// Adds a layer for the specified borehole profile based on the information given in the stringlist
	static int addLayer(std::list<std::string> fields, StationBorehole* borehole);

	/// Creates fake stratigraphies of only one layer with a thickness equal to the borehole depth
	static void createSurrogateStratigraphies(std::vector<Point*>* boreholes);

	/// Reads the specified file containing borehole stratigraphies into an vector of stringlists
	static int readStratigraphyFile(const std::string &path,
	                                std::vector<std::list<std::string> > &data);

	//long profile_type;
	//std::vector<long> _soilType;
	double _depth; // depth of the borehole
	int _date; // date when the borehole has been drilled

	/// Contains the names for all the soil layers
	std::vector<std::string> _soilName;

	/// Contains the points for the lower boundaries of all layers
	std::vector<Point*> _profilePntVec;
};
} // namespace

#endif // GEO_STATIONBOREHOLE_H
