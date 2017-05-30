/**
 * \file
 * \author Karsten Rink
 * \date   2012-08-20
 * \brief  Definition of the ImportFileTypes enumeration.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>

/**
 * \brief Types of supported import file formats.
 */
class ImportFileType
{
public:

    enum type {
        OGS = 0,
        OGS_GEO,
        OGS_STN,
        OGS_MSH,
        FEFLOW,
        GMS,
        GMSH,
        NETCDF,
        PETREL,
        POLYRASTER,
        RASTER,
        SHAPE,
        TETGEN,
        VTK
    };

    static std::string convertImportFileTypeToString(ImportFileType::type t)
    {
        if (t == ImportFileType::FEFLOW)
            return "FEFLOW";
        if (t == ImportFileType::GMS)
            return "GMS";
        if (t == ImportFileType::GMSH)
            return "GMSH";
        if (t == ImportFileType::NETCDF)
            return "NetCDF";
        if (t == ImportFileType::OGS)
            return "OGS";
        if (t == ImportFileType::OGS_GEO)
            return "OGS geometry";
        if (t == ImportFileType::OGS_STN)
            return "OGS station list";
        if (t == ImportFileType::OGS_MSH)
            return "OGS mesh";
        if (t == ImportFileType::PETREL)
            return "Petrel";
        if ((t == ImportFileType::RASTER) || (t == ImportFileType::POLYRASTER))
            return "Raster";
        if (t == ImportFileType::SHAPE)
            return "Shape";
        else if (t == ImportFileType::TETGEN)
            return "TetGen node";
        else if (t == ImportFileType::VTK)
            return "VTK";
        else
            return "";
    }

    static std::string getFileSuffixString(ImportFileType::type t)
    {
        if (t == ImportFileType::FEFLOW)
            return "FEFLOW files (*.fem)";
        if (t == ImportFileType::GMS)
            return "GMS files (*.txt *.3dm)";
        if (t == ImportFileType::GMSH)
            return "GMSH mesh files (*.msh)";
        if (t == ImportFileType::NETCDF)
            return "NetCDF files (*.nc)";
        if (t == ImportFileType::OGS)
            return "OpenGeosys files (*.gsp *.gml *.vtu *.stn);;GeoSys legacy "
                   "files (*.gli *.msh);;All files (* *.*)";
        if (t == ImportFileType::OGS_GEO)
            return "OpenGeosys files (*.gml *.gli)";
        if (t == ImportFileType::OGS_STN)
            return "OpenGeosys files (*.stn)";
        if (t == ImportFileType::OGS_MSH)
            return "OpenGeosys files (*.vtu *.msh)";
        if (t == ImportFileType::PETREL)
            return "Petrel files (*)";
        if (t == ImportFileType::RASTER)
            return "Raster files (*.asc *.grd *.bmp *.jpg *.png *.tif)";
        if (t == ImportFileType::SHAPE)
            return "ESRI Shape files (*.shp)";
        else if (t == ImportFileType::TETGEN)
            return "TetGen node files (*.node *.poly *.smesh)";
        else if (t == ImportFileType::VTK)
            return "VTK files (*.vtk *.vti *.vtr *.vts *.vtp *.vtu)";
        else
            return "All files (*.*)";
    }
};
