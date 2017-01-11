/**
 * \file
 * \author Karsten Rink
 * \date   2010-02-24
 * \brief  Definition of the VtkStationSource class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vtkPolyDataAlgorithm.h>

#include "Station.h"
#include "Applications/DataHolderLib/Color.h"
#include "VtkAlgorithmProperties.h"

/**
 * \brief VTK source object for the visualisation of station data (including boreholes)
 */
class VtkStationSource : public vtkPolyDataAlgorithm, public VtkAlgorithmProperties
{
public:
    /// Create new objects with New() because of VTKs object reference counting.
    static VtkStationSource* New();

    vtkTypeMacro(VtkStationSource,vtkPolyDataAlgorithm);

    /// Returns the colour lookup table generated for boreholes.
    /// This method should only be called after the colour lookup table has actually been build (via RequestData() or setColorLookupTable()).
    const std::map<std::string, DataHolderLib::Color>& getColorLookupTable() const
        { return _colorLookupTable; }

    /// Returns the type of observation site represented in this source object
    GeoLib::Station::StationType getType() const
        { return static_cast<GeoLib::Station*>((*_stations)[0])->type(); };

    /// Sets a predefined color lookup table for the colouring of borehole stratigraphies
//    int setColorLookupTable(const std::string &filename)
//        { return GeoLib::readColorLookupTable(_colorLookupTable, filename); }

    /// Sets the stations as a vector
    void setStations(const std::vector<GeoLib::Point*>* stations) { _stations = stations; }

    /// Prints its data on a stream.
    void PrintSelf(ostream& os, vtkIndent indent) override;

    virtual void SetUserProperty(QString name, QVariant value) override;

protected:
    VtkStationSource();

    /// Computes the polygonal data object.
    int RequestData(vtkInformation* request,
                    vtkInformationVector** inputVector,
                    vtkInformationVector* outputVector) override;

    int RequestInformation(vtkInformation* request,
                           vtkInformationVector** inputVector,
                           vtkInformationVector* outputVector) override;

    /// The stations to visualize
    const std::vector<GeoLib::Point*>* _stations;

    /// The colour table for stratigraphic data. This table is either set using the setColorLookupTable() method or is generated
    /// automatically with random colours while creating the VtkStationSource-object.
    std::map<std::string, DataHolderLib::Color> _colorLookupTable;

private:
    std::size_t GetIndexByName( std::string const& name );

    std::map<std::string, vtkIdType> _id_map;
};
