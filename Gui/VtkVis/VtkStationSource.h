/**
 * \file VtkStationSource.h
 * 24/02/2010 KR Initial implementation
 *
 */

#ifndef VTKSTATIONSOURCE_H
#define VTKSTATIONSOURCE_H

// ** INCLUDES **
#include "VtkAlgorithmProperties.h"
#include <vtkPolyDataAlgorithm.h>

#include "GEOObjects.h"

/**
 * \brief VTK source object for the visualisation of station data (including boreholes)
 */
class VtkStationSource : public vtkPolyDataAlgorithm, public VtkAlgorithmProperties
{
public:
	/// Create new objects with New() because of VTKs object reference counting.
	static VtkStationSource* New();

	vtkTypeRevisionMacro(VtkStationSource,vtkPolyDataAlgorithm);

	/// Returns the colour lookup table generated for boreholes.
	/// This method should only be called after the colour lookup table has actually been build (via RequestData() or setColorLookupTable()).
	const std::map<std::string,
	               GeoLib::Color*>& getColorLookupTable() const { return _colorLookupTable; }

	/// Sets a predefined color lookup table for the colouring of borehole stratigraphies
	int setColorLookupTable(const std::string &filename) { return readColorLookupTable(
		                                                              _colorLookupTable,
		                                                              filename); }

	/// Sets the stations as a vector
	void setStations(const std::vector<GeoLib::Point*>* stations) { _stations = stations; }

	/// Prints its data on a stream.
	void PrintSelf(ostream& os, vtkIndent indent);

	virtual void SetUserProperty(QString name, QVariant value);

protected:
	VtkStationSource();
	~VtkStationSource();

	/// Computes the polygonal data object.
	int RequestData(vtkInformation* request,
	                vtkInformationVector** inputVector,
	                vtkInformationVector* outputVector);

	int RequestInformation(vtkInformation* request,
	                       vtkInformationVector** inputVector,
	                       vtkInformationVector* outputVector);

	/// The stations to visualize
	const std::vector<GeoLib::Point*>* _stations;

	/// The colour table for stratigraphic data. This table is either set using the setColorLookupTable() method or is generated
	/// automatically with random colours while creating the VtkStationSource-object.
	std::map<std::string, GeoLib::Color*> _colorLookupTable;

private:
	size_t GetIndexByName( std::string name );

	std::map<std::string, vtkIdType> _id_map;
};

#endif // VTKSTATIONSOURCE_H
