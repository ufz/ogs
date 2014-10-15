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

std::vector<std::string> split(const std::string &str, char delim)
{
    std::vector<std::string> res;
    std::size_t current = 0, found;
    while ((found = str.find_first_of(delim, current)) != std::string::npos) {
        std::string sub(str, current, found - current);
        BaseLib::trim(sub,'"');
        res.push_back(sub);
        current = found + 1;
    }
    std::string sub(str, current, str.size() - current);
    BaseLib::trim(sub,'"');
    res.push_back(sub);
    return res;
}

unsigned findInSplit(const std::vector<std::string> splitted, const std::string key)
{
    unsigned id = 0;
    for (auto s : splitted) {
        if (s.compare(key)==0) return id;
        id++;
    }
    return -1;
}

void readParaViewCSV(const std::string &filename, RasterData &raster)
{
    std::ifstream is(filename.c_str());
    if (!is.is_open()) {
        std::cout << "***ERROR: cannot open a file " << filename << std::endl;
        exit(1);
    }

    raster.ncols = 1780;
    raster.nrows = 919;
    raster.cellsize = 0.05;
    raster.xllcorner = .0;
    raster.yllcorner = .0;
    raster.NODATA_Value = -9999;
    raster.values.resize(raster.ncols*raster.nrows);
    std::cout << "-> expecting " << raster.values.size() << " points" << std::endl;

    std::string buffer;
    std::size_t n_lines = 0;

    // header
    std::getline(is, buffer);
    const auto columns = split(buffer, ',');
    n_lines++;
    std::cout << "-> found " << columns.size() << " columns in this CSV" << std::endl;
    for (auto key : columns) std::cout << "\t" << key << std::endl;
    const unsigned colId_x = findInSplit(columns, "Points:0");
    const unsigned colId_y = findInSplit(columns, "Points:1");
    const unsigned colId_b = findInSplit(columns, "FRACTURE_GAP_NORMAL");

    while (is && std::getline(is, buffer)) {
        const auto columns = split(buffer, ',');
        n_lines++;
        double x = BaseLib::str2number<double>(columns[colId_x]);
        double y = BaseLib::str2number<double>(columns[colId_y]);
        double b = BaseLib::str2number<double>(columns[colId_b]);
        size_t col = x / raster.cellsize;
        size_t row = y / raster.cellsize;
        size_t cellid = col + row*raster.ncols;
        raster.values[cellid] = b;
    }
    std::cout << "-> the CSV file had " << n_lines << " lines" << std::endl;

}

void outputAscii(const RasterData &raster, const std::string &filename)
{
    std::ofstream os(filename.c_str());
    if (!os.is_open()) {
        std::cout << "***ERROR: cannot open a file " << filename << std::endl;
        exit(1);
    }

    os << " ncols " << raster.ncols << "\n";
    os << " nrows " << raster.nrows << "\n";
    os << " xllcorner " << raster.xllcorner << "\n";
    os << " yllcorner " << raster.yllcorner << "\n";
    os << " cellsize " << raster.cellsize << "\n";
    os << " NODATA_Value " << raster.NODATA_Value << "\n";
    os.precision(std::numeric_limits<double>::digits10);
    for (std::size_t i=0; i<raster.nrows; i++) {
        for (std::size_t j=0; j<raster.ncols; j++) {
            os << raster.values[(raster.nrows-i-1)*raster.ncols+j] << " "; // inverse order in row
        }
        os << "\n";
    }

    os.close();
}

void checkRaster(RasterData &raster)
{
    size_t n_zero=0;
    size_t n_negative = 0;
    for (size_t i=0; i<raster.values.size(); i++) {
        if (raster.values[i]==0) n_zero++;
        else if (raster.values[i]<0) {
            raster.values[i] = 1e-6; // 1nm
            n_negative++;
        }
    }
    std::cout << "-> " << n_zero << " cells have zero" << std::endl;
    std::cout << "-> " << n_negative << " cells have negative values" << std::endl;
}

int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cout << "USAGE:" << std::endl;
        std::cout << "\t" << argv[0] << " (input mesh file) (output ASCII file)" << std::endl;
        return 0;
    }
    std::string in_filename(argv[1]);
    std::string out_filename(argv[2]);

    RasterData data;
    std::cout << "-> reading a ParaView CSV file..." << std::endl;
    readParaViewCSV(in_filename, data);
    std::cout << "-> checking the data..." << std::endl;
    checkRaster(data);
    std::cout << "-> writing a raster ASCII file..." << std::endl;
    outputAscii(data, out_filename);

	return 0;
}
