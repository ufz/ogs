/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef FEFLOWMESHINTERFACE_H_
#define FEFLOWMESHINTERFACE_H_

#include <iosfwd>
#include <string>
#include <vector>

namespace GeoLib
{
class Point;
class Polyline;
}

namespace MeshLib
{
class Mesh;
class Element;
class Node;
enum class MeshElemType;
}

namespace FileIO
{
/**
 * Read FEFLOW model files (*.fem) into OGS data structure. Currently this class
 * supports only import of mesh data and some geometry given in Supermesh section.
 */
class FEFLOWMeshInterface
{
public:
    /**
     * read a FEFLOW Model file (*.fem) in ASCII format (Version 5.4)
     *
     * This function reads mesh data in addition to geometry data given in
     * Supermesh.
     *
     * @param filename  FEFLOW file name
     * @return a pointer to a created OGS mesh
     */
    MeshLib::Mesh* readFEFLOWFile(const std::string& filename);

private:
    // CLASS
    struct FEM_CLASS
    {
        unsigned problem_class = 0;
        unsigned time_mode = 0;
        unsigned orientation = 0;
        unsigned dimension = 0;
        unsigned n_layers3d = 0;
        unsigned saturation_flag = 0;
        unsigned save_fsize_rreal = 0;
        unsigned save_fsize_creal = 0;
    };

    // DIMENSION
    struct FEM_DIM
    {
        std::size_t n_nodes = 0;
        std::size_t n_elements = 0;
        std::size_t obs = 0;
        std::size_t np_cor = 0;
        unsigned n_nodes_of_element = 0;
        unsigned n_steps = 0;
        unsigned icrank = 0;
        unsigned upwind = 0;
        unsigned optim = 0;
        unsigned aquifer_type = 0;
        unsigned nwca = 0;
        unsigned adaptive_mesh = 0;
        unsigned sp_fem_pcs_id = 0;
        unsigned sorption_type = 0;
        unsigned reaction_type = 0;
        unsigned dispersion_type = 0;
    };

    /// Read element type and node indices according to the element type.
    MeshLib::Element* readElement(std::string const& line,
                                  std::vector<MeshLib::Node*> const& nodes);

    /// read node indices and create a mesh element
    MeshLib::Element* readElement(const FEM_DIM& fem_dim,
                                  const MeshLib::MeshElemType elem_type,
                                  const std::string& line,
                                  const std::vector<MeshLib::Node*>& nodes);

    /// read node coordinates
    void readNodeCoordinates(std::ifstream& in,
                             const FEM_CLASS& fem_class,
                             const FEM_DIM& fem_dim,
                             std::vector<MeshLib::Node*>& nodes);

    /// read elevation data
    void readElevation(std::ifstream& in,
                       const FEM_CLASS& fem_class,
                       const FEM_DIM& fem_dim,
                       std::vector<MeshLib::Node*>& vec_nodes);

    //// parse node lists
    std::vector<std::size_t> getIndexList(const std::string& str_ranges);

    /// parse ELEMENTALSETS
    void readELEMENTALSETS(
        std::ifstream& in,
        std::vector<std::vector<std::size_t>>& vec_elementsets);

    void setMaterialIDs(
        FEM_CLASS const& fem_class,
        FEM_DIM const& fem_dim,
        std::vector<GeoLib::Polyline*>* const& lines,
        std::vector<std::vector<std::size_t>> const& vec_elementsets,
        std::vector<MeshLib::Element*> const& vec_elements,
        std::vector<int>& material_ids);
};
}  // FileIO

#endif /* FEFLOWMESHINTERFACE_H_ */
