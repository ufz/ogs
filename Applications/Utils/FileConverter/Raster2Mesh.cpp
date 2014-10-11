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

#if 0
void outputTriMeshSubdivide(const RasterData &raster, const std::string &filename)
{
    const std::size_t n_org_cells = raster.ncols * raster.nrows;
    const std::size_t n_div = 2;
    const std::size_t n_new_eles = n_org_cells*n_div*n_div*2;
    const std::size_t n_new_nodes = (raster.ncols*n_div+1) * (raster.nrows*n_div+1);

    std::ofstream os(filename.c_str());
    os.precision(std::numeric_limits<double>::digits10);
    os << "#FEM_MSH\n";
    os << " $PCS_TYPE\n";
    os << "  NO_PCS\n";
    os << " $NODES\n";
    os << "  " << n_new_nodes << "\n";
    std::size_t node_id = 0;
    const double dL = raster.cellsize / n_div;
    double z = .0;
    std::vector<std::size_t> cells(4);
    for (std::size_t i=0; i<raster.nrows*n_div+1; i++) {
        for (std::size_t j=0; j<raster.ncols*n_div+1; j++) {
            int type = 0;
            const std::size_t cell_id = (i/n_div)*raster.ncols + (j/n_div);
            if ((i+1)%2==0 && (j+1)%2==0) {
                type=0; // cell center
            }  else {
                if ((i==0 || i==raster.nrows*n_div) && (j==0 || j==raster.ncols*n_div) ) {
                    type=0; // corner points, same value as cell center
                } else if (i==0 || i==raster.nrows*n_div) {
                    if ((j+1)%2==0) {
                        type = 0; // lower or upper edge, in cell
                    } else {
                        type= 1; // lower or lower edge, value should be interpolated
                        cells[0] = (i/n_div)*raster.ncols + j/n_div-1;
                        cells[1] = (i/n_div)*raster.ncols + j/n_div;
                    }
                } else if (j==0 || j==raster.ncols*n_div) {
                    if ((i+1)%2==0) {
                        type = 0; // left or right edge, in cell
                    } else {
                        type= 1; // left or right edge, value should be interpolated
                        cells[0] = (i/n_div-1)*raster.ncols + j/n_div;
                        cells[1] = (i/n_div)*raster.ncols + j/n_div;
                    }
                } else if (i%n_div==0 && j%n_div==0) {
                    // 4 neighbors
                    type = 2;
                } else {
                    // 2 neighbors
                    type = 1;
                }

            }
            //
            switch (type)
            {
            case 0:
                {
                    cells[0] = (i/n_div)*raster.nrows + (j/n_div);
                    z = raster.values[cells[0]];
                }
                break;
            case 1:
                {

                }
                break;
            }
            os << "  " << node_id++ << j*dL+raster.xllcorner << " " <<  i*dL+raster.yllcorner << " " << z << "\n";
        }
    }
    os << " $ELEMENTS\n";
    os << "  " << n_new_eles << "\n";
    std::size_t ele_id = 0;
    std::size_t n_row_nodes = raster.nrows*n_div + 1;
    for (std::size_t i=0; i<raster.nrows*n_div; i++) {
        for (std::size_t j=0; j<raster.ncols*n_div; j++) {
            os << "  " << ele_id++ << " 0 tri "  << i*n_row_nodes + j << " " <<  i*n_row_nodes + j +1 << " " << (i+1)*n_row_nodes + j + 1 << "\n";
            os << "  " << ele_id++ << " 0 tri "  << i*n_row_nodes + j << " " <<  (i+1)*n_row_nodes + j + 1 << " " << (i+1)*n_row_nodes + j << "\n";
        }
    }
    os << "#STOP\n";
    os.close();
}
#endif

void outputTriMesh(const RasterData &raster, const std::string &filename)
{
    const std::size_t n_org_cells = raster.ncols * raster.nrows;
    const std::size_t n_new_nodes = n_org_cells;
    const std::size_t n_new_eles = (raster.ncols-1)*(raster.nrows-1)*2;

    std::ofstream os(filename.c_str());
    os.precision(std::numeric_limits<double>::digits10);
    os << "#FEM_MSH\n";
    os << " $PCS_TYPE\n";
    os << "  NO_PCS\n";
    os << " $NODES\n";
    os << "  " << n_new_nodes << "\n";
    std::size_t node_id = 0;
    const double dL = raster.cellsize;
    const double offset_x = raster.xllcorner+dL*0.5;
    const double offset_y = raster.yllcorner+dL*0.5;
    for (std::size_t i=0; i<raster.nrows; i++) {
        for (std::size_t j=0; j<raster.ncols; j++) {
        	// output x,y from lower left to upper right
        	// Remark: values are stored from upper left to bottom right
            os << "  " << node_id++ << " " << j*dL+offset_x << " " <<  i*dL+offset_y << " " << raster.values[(raster.nrows-i-1)*raster.ncols+j] << "\n";
        }
    }
    os << " $ELEMENTS\n";
    os << "  " << n_new_eles << "\n";
    std::size_t ele_id = 0;
    for (std::size_t i=0; i<raster.nrows-1; i++) {
        for (std::size_t j=0; j<raster.ncols-1; j++) {
            os << "  " << ele_id++ << " 0 tri "  << i*raster.ncols + j << " " <<  i*raster.ncols + j +1 << " " << (i+1)*raster.ncols + j + 1 << "\n";
            os << "  " << ele_id++ << " 0 tri "  << i*raster.ncols + j << " " <<  (i+1)*raster.ncols + j + 1 << " " << (i+1)*raster.ncols + j << "\n";
        }
    }
    os << "#STOP\n";
    os.close();
}

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

int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cout << "USAGE:" << std::endl;
        std::cout << "\t" << argv[0] << " (input ASCII file) (output mesh file)" << std::endl;
        return 0;
    }
    std::string in_filename(argv[1]);
    std::string out_filename(argv[2]);

    RasterData data;
    std::cout << "-> reading a raster file..." << std::endl;
    readRaster(in_filename, data);
    std::cout << "-> creating a triangle mesh file..." << std::endl;
    outputTriMesh(data, out_filename);

	return 0;
}
