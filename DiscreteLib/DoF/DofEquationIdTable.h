/**
 * \file   DofEquationIdTable.h
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef DOFEQUATIONIDTABLE_H_
#define DOFEQUATIONIDTABLE_H_

#include <vector>
#include <map>

#include "BaseLib/CodingTools.h"

#include "DofNumberingType.h"
#include "IEquationIdStorage.h"
#include "DoF.h"

namespace DiscreteLib
{

/**
 * \brief Mapping table between DoFs and equation index
 *
 * This class contains a table which relates degree of freedoms (DoFs) and
 * indexes in a system of equations (Equation ID).
 *
 * Equation IDs can be numbered either by variable id or by Point ID of DoFs.
 * DofEquationIdTable class can have different numbering types for global and
 * local view.
 *
 * Example Usage:
 * \code
 *   std::size_t mesh_id = 0;
 *   DofEquationIdTable dofMap;
 *   std::size_t var_id = dofMap.addVariableDoFs(mesh_id, 0, 10); // add 10 DoFs from point 0
 *   dofMap.setNumberingType(DofNumberingType::BY_POINT);
 *   dofMap.construct();
 *   std::size_t pt0_eqsId = dofMap.mapEqsID(DoF(var_id, mesh_id, 0));
 * \endcode
 */
class DofEquationIdTable
{
public:

    ///
    DofEquationIdTable()
    : _total_dofs(0), _total_dofs_without_ghosts(0), _is_seq_address(false),
      _global_numbering_type(DofNumberingType::BY_POINT), _local_numbering_type(DofNumberingType::BY_VARIABLE)
    {};

    ///
    virtual ~DofEquationIdTable();

    //---------------------------------------------------------------------------------------------
    // setup DoFs
    //---------------------------------------------------------------------------------------------
    /**
     * add a variable and related DoFs with random point IDs
     *
     * @param[in] mesh_id       Mesh ID
     * @param[in] vec_pt_id     Vector of discrete point id
     * @return Variable ID
     */
    std::size_t addVariableDoFs(std::size_t mesh_id, const std::vector<std::size_t> &vec_pt_id);

    /**
     * add a variable and related DoFs with sequential point IDs
     *
     * @param mesh_id           Mesh ID
     * @param pt_id_begin   the beginning id of discrete points
     * @param pt_size      the number of discrete points
     * @return Variable ID
     */
    std::size_t addVariableDoFs(std::size_t mesh_id, std::size_t pt_id_begin, std::size_t pt_size);

    //---------------------------------------------------------------------------------------------
    // construct the table
    //---------------------------------------------------------------------------------------------
    /**
     * set numbering type
     * @param numbering
     */
    void setNumberingType(DofNumberingType::type numbering) {_global_numbering_type = numbering;};

    /**
     * get numbering type
     * @return
     */
    DofNumberingType::type getNumberingType() const {return _global_numbering_type;};

    /**
     * set local numbering type
     * @param numbering
     */
    void setLocalNumberingType(DofNumberingType::type numbering) {_local_numbering_type = numbering;};

    /**
     * get local numbering type
     * @return
     */
    DofNumberingType::type getLocalNumberingType() const {return _local_numbering_type;};

    /**
     * construct the DoF mapping table
     *
     * The followings will be done by this function:
     * - numbering DoFs
     * - count the total number of DoFs with/without ghost point
     *
     * @param[in] offset
     */
    virtual void construct(std::size_t offset=0);

    //---------------------------------------------------------------------------------------------
    // get summary of the table
    //---------------------------------------------------------------------------------------------
    //# Mesh
    /// get the number of meshes
    std::size_t getNumberOfMeshes() const {return _map_msh2var.size();};

    //# Variables
    /// get the number of registered variables
    std::size_t getNumberOfVariables() const { return _map_var2dof.size(); }

    /**
     * get the number of variables associated with the given mesh
     * @param mesh_id   Mesh ID
     * @return
     */
    std::size_t getNumberOfVariables(std::size_t mesh_id) const
    {
        auto itr = _map_msh2var.find(mesh_id);
        return (itr!=_map_msh2var.end()) ? itr->second.size() : 0;
    };

    //---------------------------------------------------------------------------------------------
    // Active and inactive DoFs
    //---------------------------------------------------------------------------------------------
    //# DoFs
    /**
     * activate or deactivate DoF at the given point
     *
     * @param dof       DoF
     * @param active    Active or not
     */
    void activateDoF(const DoF &dof, bool active);

    /**
     * is this DoF active?
     * @param dof       DoF
     * @return
     */
    bool isActiveDoF(const DoF &dof) const;

    /// get the total number of active DoFs
    std::size_t getTotalNumberOfActiveDoFs() const { return _total_dofs; }

    /// get the total number of active DoFs excluding ghost points
    std::size_t getTotalNumberOfActiveDoFsWithoutGhost() const {return _total_dofs_without_ghosts;};

    //---------------------------------------------------------------------------------------------
    // Ghost points
    //---------------------------------------------------------------------------------------------
    //# Ghost 
    /**
     * set ghost points for overlapped regions appeared with domain-decomposition methods
     * @param mesh_id       Mesh ID
     * @param list_pt_id    Vector of ghost point id
     */
    void setGhostPoints(std::size_t mesh_id, const std::vector<std::size_t> &list_pt_id);

    /**
     * is this point a ghost?
     * @param mesh_id   Mesh ID
     * @param pt_id     Discrete point ID
     * @return
     */
    bool isGhostPoint(std::size_t mesh_id, std::size_t pt_id) const;

    /**
     * get the number of ghost points
     * @param mesh_id   Mesh ID
     * @return
     */
    std::size_t getNumberOfGhostPoints(std::size_t mesh_id) const
    {
        auto itr = _ghost_pt.find(mesh_id);
        return (itr!=_ghost_pt.end()) ? itr->second.size() : 0;
    }

    /// get the total number of ghost points
    std::size_t getTotalNumberOfGhostPoints() const;

    //---------------------------------------------------------------------------------------------
    // mapping between DoF and Equation ID
    //---------------------------------------------------------------------------------------------
    /**
     * return an equation id for the given DoF
     *
     * @param dof       DoF
     * @return equation id or npos if no match
     */
    std::size_t mapEqsID(const DoF &dof) const;

    /**
     * return a list of equation id for the given variable, mesh and points
     *
     * @param[in] var_id        Variable ID
     * @param[in] mesh_id       Mesh ID
     * @param[in] vec_pt_id     Vector of point ID
     * @param[out] vec_eqs_id   Vector of equation ID
     */
    void mapEqsID(  std::size_t var_id,
                    std::size_t mesh_id,
                    const std::vector<std::size_t> &vec_pt_id,
                    std::vector<std::size_t> &vec_eqs_id) const;

    /**
     * return a list of equation id for all variables, the given mesh and points
     *
     * @param[in] mesh_id       Mesh ID
     * @param[in] vec_pt_id     Vector of point ID
     * @param[out] vec_eqs_id   Vector of equation ID
     *  The length of the returned vector is the number of points multiplied by the
     *  number of variables. Equation IDs are sorted by the given local numbering
     *  type.
     */
    void mapEqsID(  std::size_t mesh_id,
                    const std::vector<std::size_t> &vec_pt_id,
                    std::vector<std::size_t> &vec_eqs_id) const;

    /**
     * return a list of equation id for all variables, the given mesh and points
     *
     * @param[in] mesh_id                   Mesh ID
     * @param[in] pt_id                     Vector of point ID
     * @param[out] vec_eqs_id                    Vector of equation ID
     *  The length of the returned vector is the number of points multiplied by the
     *  number of variables. Equation IDs are sorted by the given local numbering
     *  type. A created list may contain entries with -1 if no valid address is found
     *  for corresponding key (var id, point id).
     * @param[out] vec_eqs_id_without_ghost      Vector of equation ID without ghost point
     *  The vector may have -1 also for ghost points.
     */
    void mapEqsID(  std::size_t mesh_id,
                    const std::vector<std::size_t> &list_pt_id,
                    std::vector<std::size_t> &vec_eqs_id,
                    std::vector<std::size_t> &vec_eqs_id_without_ghost ) const;
    
    /**
     * return a list of equation id for all variables
     *
     * Length of a created list is sum of the number of active points for each variable.
     * The second list may have -1 for ghost points.
     *
     * @param mesh_id
     * @param pt_id
     * @param eqs_id
     * @param eqs_id_without_ghost
     */
    void mapEqsIDreduced(   std::size_t mesh_id,
                            const std::vector<std::size_t> &list_pt_id,
                            std::vector<std::size_t> &list_eqs_id,
                            std::vector<std::size_t> &list_eqs_id_without_ghost
                            ) const;

    /**
     * return a DoF from the given equation ID
     *
     * @param eqs_id    Equation ID
     * @return DoF
     */
    DoF mapDoF(std::size_t eqs_id) const;

    //---------------------------------------------------------------------------------------------
    // some helper functions
    //---------------------------------------------------------------------------------------------
    /**
     * create a local mapping table
     *
     * This function is used for conversion between local and global linear systems.
     *
     * @param mesh_id       Mesh ID
     * @param vec_pt_id     Vector of point ID
     * @param local_table   Local mapping table
     *  A created table contains mapping relations between the given point IDs but
     *  local equation IDs (local means only DoFs with the given points are considered).
     */
    void createLocalMappingTable(   std::size_t mesh_id,
                                    const std::vector<std::size_t> &vec_pt_id,
                                    DofEquationIdTable &local_table
                                    ) const;

    /// print debug info
    void printout() const;

protected:
    /**
     * add a variable with the given mapping
     *
     * @param mesh_id
     * @param address
     * @return
     */
    std::size_t addVariable(std::size_t mesh_id, IEquationIdStorage* address);

    /**
     *
     * @param var_id
     * @param mesh_id
     * @return
     */
    const IEquationIdStorage* getEquationIdStorage(std::size_t var_id, std::size_t mesh_id) const;

    /**
     *
     * @param var_id
     * @param mesh_id
     * @return
     */
    IEquationIdStorage* getEquationIdStorage(std::size_t var_id, std::size_t mesh_id) { return _map_var2dof[var_id][mesh_id]; }

private:
    DISALLOW_COPY_AND_ASSIGN(DofEquationIdTable);


private:
    typedef std::size_t VariableID;
    typedef std::size_t MeshID;
    typedef std::size_t PointID;
    std::vector<std::map<MeshID, IEquationIdStorage*> > _map_var2dof; // var -> (mesh, address)
    std::map<MeshID, std::vector<VariableID> > _map_msh2var; // mesh -> (var)
    std::map<MeshID, std::vector<PointID> > _ghost_pt; // mesh -> ghost points

    std::size_t _total_dofs;
    std::size_t _total_dofs_without_ghosts;
    bool _is_seq_address;

    DofNumberingType::type _global_numbering_type;
    DofNumberingType::type _local_numbering_type;
};

} //end

#endif //DOFEQUATIONIDTABLE_H_
