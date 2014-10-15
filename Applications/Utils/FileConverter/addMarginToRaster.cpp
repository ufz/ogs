//============================================================================
// Name        : fracture_map.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>

#include "BaseLib/StringTools.h"

struct RasterData
{
    std::size_t ncols;
    std::size_t nrows;
    double xllcorner;
    double yllcorner;
    double cellsize;
    double NODATA_Value;
    std::vector<double> values;
};

void readRaster(const std::string &filename, RasterData &raster)
{
    std::ifstream is(filename.c_str());
    std::string dummy;
    is >> dummy >> raster.ncols;
    is >> dummy >> raster.nrows;
    is >> dummy >> raster.xllcorner;
    is >> dummy >> raster.yllcorner;
    is >> dummy >> raster.cellsize;
    is >> dummy >> raster.NODATA_Value;
    const std::size_t n_data = raster.ncols * raster.nrows;
    raster.values.resize(n_data);
    for (std::size_t i=0; i<n_data; i++)
        is >> raster.values[i];

    is.close();
}

void writeRaster(const RasterData &raster, const std::string &filename)
{
    std::ofstream os(filename.c_str());
    if (!os.is_open()) {
        std::cout << "***ERROR: cannot open a file " << filename << std::endl;
        exit(1);
    }

    os << "ncols " << raster.ncols << "\n";
    os << "nrows " << raster.nrows << "\n";
    os << "xllcorner " << raster.xllcorner << "\n";
    os << "yllcorner " << raster.yllcorner << "\n";
    os << "cellsize " << raster.cellsize << "\n";
    os << "NODATA_Value " << raster.NODATA_Value << "\n";
    os.precision(std::numeric_limits<double>::digits10);
    for (std::size_t i=0; i<raster.nrows; i++) {
        for (std::size_t j=0; j<raster.ncols; j++) {
            os << raster.values[i*raster.ncols+j] << " "; // inverse order in row
//            os << raster.values[(raster.nrows-i-1)*raster.ncols+j] << " "; // inverse order in row
        }
        os << "\n";
    }

    os.close();
}

void addMerginToRaster(RasterData &raster, const double dX)
{
    const size_t dColumns = dX/raster.cellsize;
    raster.ncols += 2*dColumns;
    raster.xllcorner = -dX;
    std::vector<double> oldvalues(raster.values);

    raster.values.resize(raster.nrows*raster.ncols);
    for (std::size_t i=0; i<raster.nrows; i++) {
        for (std::size_t j=0; j<raster.ncols; j++) {
            if (j>dColumns-1 && j<raster.ncols-dColumns)
                raster.values[i*raster.ncols+j] = oldvalues[i*(raster.ncols-2*dColumns)+j-dColumns];
            else
                raster.values[i*raster.ncols+j] = 10.0; // 10mm
        }
    }
}

int main(int argc, char* argv[])
{
    if (argc < 4) {
        std::cout << "USAGE:" << std::endl;
        std::cout << "\t" << argv[0] << " (input ASCII file) (output ASCII file) (dx)" << std::endl;
        return 0;
    }
    std::string in_filename(argv[1]);
    std::string out_filename(argv[2]);
    const double dx = BaseLib::str2number<double>(argv[3]);

    RasterData data;
    std::cout << "-> reading a raster file..." << std::endl;
    readRaster(in_filename, data);
    std::cout << "-> checking the data..." << std::endl;
    addMerginToRaster(data, dx);
    std::cout << "-> writing a raster ASCII file..." << std::endl;
    writeRaster(data, out_filename);

	return 0;
}
