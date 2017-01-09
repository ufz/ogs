/**
 * \file
 * \author Karsten Rink
 * \date   2010-04-23
 * \brief  Definition of the VtkColorLookupTable class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VTKCOLORLOOKUPTABLE_H
#define VTKCOLORLOOKUPTABLE_H

// ** INCLUDES **
#include <cassert>
#include <map>
#include <vector>
#include <vtkLookupTable.h>

#include "Applications/DataHolderLib/ColorLookupTable.h"

/**
 * \brief Calculates and stores a colour lookup table.
 *
 * Based on a start colour and an end colour, RGB-values are interpolated and stored in vector of GeoLib::Color. If no
 * colours are set, default values are used for start (blue) and end (red). The number of entries of the colour table can
 * be set in the constructor, the default value is 256. If additional colours are inserted into the table using setColor()
 * the interpolation will be calculated iteratively between set colour values. Interpolation can be linear (default) or
 * exponential. Based on the set range of values, colour values can be retrieved using getColor().
 */
class VtkColorLookupTable : public vtkLookupTable
{
public:
    static const int DEFAULTMINVALUE = -9999;
    static const int DEFAULTMAXVALUE =  9999;

    /// \brief Create new objects with New() because of VTKs object reference counting.
    static VtkColorLookupTable* New();

    vtkTypeMacro(VtkColorLookupTable,vtkLookupTable);

    /// \brief Builds the colour table based on the previously set parameters.
    /// This method should only be called after all options have been set.
    void Build();

    /* \brief Sets the given colour as a constant in the colour lookup table.
     * The colour will subsequently be considered in the interpolation process when the lookup table is built. Note that pos is only a
     * relative position, i.e. pos in (0,1). The actual position of that colour in the table is dependent on the number of entries set
     * the SetRange() method.
     */
    void setColor(double pos, DataHolderLib::Color const& color);

    /* \brief Returns the colour at the given index from the colour lookup table.
     * The colour will be interpolated from the colour-dictionary entries before and after this index.
     * Make sure that Build() has been called before using this method.
     */
    void getColor(vtkIdType indx, unsigned char rgba[4]) const;

    /// Returns the type of interpolation used.
    DataHolderLib::LUTType getInterpolationType() const { return _type; }

    /// Sets the type of interpolation.
    void setInterpolationType(DataHolderLib::LUTType type) { _type = type; }

    /// Imports settings of OGS lookup table class
    void setLookupTable(DataHolderLib::ColorLookupTable const& lut);

    /// Exports a color table to a file.
    void writeToFile(const std::string &filename);

    /// Set a value within the LUT
    void SetTableValueRGBA(vtkIdType idx, unsigned char rgba[4]);

    /// Get a value from the LUT
    void GetTableValue(vtkIdType idx, unsigned char rgba[4]);

protected:
    /// Constructor
    VtkColorLookupTable();

    /// Destructor
    ~VtkColorLookupTable();

private:
    /// Interpolates values linearly.
    unsigned char linInterpolation(unsigned char a, unsigned char b, double p) const;

    /// Interpolates values exponentially. gamma should roughly be in [0,4), for gamma=1 interpolation is linear.
    unsigned char expInterpolation(unsigned char a, unsigned char b, double gamma, double p) const;

    std::map<double, unsigned char*> _dict;
    DataHolderLib::LUTType _type;
};

#endif // VTKCOLORLOOKUPTABLE_H
