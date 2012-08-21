/**
 * \file GMSInterface.h
 * 08/06/2010 KR Initial implementation
 *
 */

#ifndef GMSINTERFACE_H_
#define GMSINTERFACE_H_

#include <list>
#include <vector>

#include "Point.h"

namespace GeoLib {
	class Station;
	class StationBorehole;
}

namespace MeshLib {
	class Mesh;
}

/**
 * \brief Manages the import and export of Aquaveo GMS files into and out of GEOLIB.
 */
class GMSInterface
{
public:
	/// Exports borehole data from all boreholes in a list to a file in GMS-format. (Note: there are some hardcoded tmp-files in the method that you might need to change!)
	static void writeBoreholesToGMS(const std::vector<GeoLib::Point*>* stations,
	                                const std::string &filename);

	/// Imports borehole data from a file in GMS-format.
	static int readBoreholesFromGMS(std::vector<GeoLib::Point*>* boreholes,
	                                const std::string &filename);

	/// Exports borehole data from one borehole to a file in GMS-format.
	static int writeBoreholeToGMS(const GeoLib::StationBorehole* station,
	                              const std::string &filename,
	                              std::vector<std::string> &soilID);

	/// Writes a file that assigns each soilID-index in the GMS export file a name.
	static int writeSoilIDTable(const std::vector<std::string> &soilID,
	                            const std::string &filename);

	/// Reads a GMS *.3dm file and converts it to an CFEMesh.
	static MeshLib::Mesh* readGMS3DMMesh(std::string file_name);

private:
	/**
	 * \brief Reads SoilIDs for Borehole export from an external file
	 *
	 * The method expects a file with the name of one stratigraphic layer at each line. These layers are assigned
	 * ascending IDs, i.e. the first name gets index 0, the second line gets index 1, etc.
	 * \return An array with the names of the stratigraphic layers in which the index for each string equals its ID.
	 */
	static std::vector<std::string> readSoilIDfromFile(const std::string &filename);

	/// Finds the ID assigned to soilName or creates a new one ( this method is called from writeBoreholeToGMS() )
	static size_t getSoilID(std::vector<std::string> &soilID, std::string &soilName);
};

#endif /* GMSINTERFACE_H_ */
