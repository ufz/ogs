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
#include <cassert>
#include <algorithm>

#include "BaseLib/CodingTools.h"

#include "DofNumberingType.h"
#include "IEquationIdStorage.h"


namespace DiscreteLib
{

/**
 * \brief Mapping table between DoFs and equation index
 *
 * This class contains a table which relates DoFs and index in a system of equations.
 * A DoF is implicitly defined as a unique combination of the followings 
 * - variable type
 * - mesh id
 * - discrete point id in the mesh (e.g. node id, element id)
 * 
 */
class DofEquationIdTable
{
public:

    ///
    DofEquationIdTable()
    : _total_dofs(0), _total_dofs_without_ghosts(0), _is_seq_address(false),
      _numbering_type(DofNumberingType::BY_POINT), _local_numbering_type(DofNumberingType::BY_VARIABLE)
    {};

    ///
    virtual ~DofEquationIdTable();

    //---------------------------------------------------------------------------------------------
    // setup DoFs
    //---------------------------------------------------------------------------------------------
    /**
     * add DoFs of a variable
     * @param mesh_id           Mesh ID
     * @param list_dof_pt_id    Vector of discrete point id
     * @return Variable ID
     */
    std::size_t addVariableDoFs(std::size_t mesh_id, const std::vector<std::size_t> &list_dof_pt_id);

    /**
     * add DoFs of a variable
     * @param mesh_id           Mesh ID
     * @param dof_pt_id_begin   the beginning id of discrete points
     * @param dof_pt_count      the number of discrete points
     * @return Variable ID
     */
    std::size_t addVariableDoFs(std::size_t mesh_id, std::size_t dof_pt_id_begin, std::size_t dof_pt_count);

    /**
     * add DoFs with multiple variables but using the same discrete points
     * @param n_var             The number of variables
     * @param mesh_id           Mesh ID
     * @param dof_pt_id_begin   Starting point ID
     * @param dof_pt_count      The total number of points
     */
    void addVariableDoFs(std::size_t n_var, std::size_t mesh_id, std::size_t dof_pt_id_begin, std::size_t dof_pt_count);


    //---------------------------------------------------------------------------------------------
    // construct the table
    //---------------------------------------------------------------------------------------------
    /**
     * set numbering type
     * @param numbering
     */
    void setNumberingType(DofNumberingType::type numbering) {_numbering_type = numbering;};

    /**
     * get numbering type
     * @return
     */
    DofNumberingType::type getNumberingType() const {return _numbering_type;};

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
     * numbering DoFs
     * @param offset
     */
    virtual void construct(long offset=0);

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

    //# DoFs
    /**
     * deactivate DoFs
     * @param var_id    Variable ID
     * @param mesh_id   Mesh ID
     * @param pt_id     Point ID to be deactivated
     */
    void deactivateDoFs(std::size_t var_id, std::size_t mesh_id, std::size_t pt_id);
    /**
     * is this DoF active?
     * @param var_id    Variable ID
     * @param mesh_id   Mesh ID
     * @param pt_id     Point ID
     * @return
     */
    bool isActiveDoF(std::size_t var_id, std::size_t mesh_id, std::size_t pt_id) const;
    /// get the total number of active DoFs
    std::size_t getTotalNumberOfActiveDoFs() const { return _total_dofs; }
    /// get the total number of active DoFs excluding ghost points
    std::size_t getTotalNumberOfActiveDoFsWithoutGhost() const {return _total_dofs_without_ghosts;};


    //---------------------------------------------------------------------------------------------
    // equation address mapping
    //---------------------------------------------------------------------------------------------
    /**
     * return equation id
     *
     * @param var_id    Variable ID
     * @param mesh_id   Mesh ID
     * @param pt_id     Point ID
     * @return equation id or npos if no match
     */
    std::size_t mapEqsID(std::size_t var_id, std::size_t mesh_id, std::size_t pt_id) const;

    /**
     * return a list of equation id
     *
     * @param var_id        Variable ID
     * @param mesh_id       Mesh ID
     * @param pt_id         Vector of point ID
     * @return vector of equation id
     */
    void mapEqsID(std::size_t var_id, std::size_t mesh_id, const std::vector<std::size_t> &pt_id, std::vector<std::size_t> &eqs_id) const;

    /**
     * return a list of equation id for all variables
     *
     * @param mesh_id   Mesh ID
     * @param pt_id     Vector of point ID
     * @param eqs_id    Vector of equation ID
     */
    void mapEqsID(std::size_t mesh_id, const std::vector<std::size_t> &pt_id, std::vector<std::size_t> &eqs_id) const;

    /**
     * return a list of equation id for all variables
     *
     * Length of a created list is the number of points * the number of variables.
     * a created list may contain entries with -1 if no valid address is found
     * for corresponding key (var id, point id). In addition to that, the second
     * may have -1 for ghost points.
     *
     * @param mesh_id                   Mesh ID
     * @param pt_id                     Vector of point ID
     * @param eqs_id                    Vector of equation ID
     * @param eqs_id_without_ghost      Vector of equation ID without ghost point
     */
    void mapEqsID(
            std::size_t mesh_id, const std::vector<std::size_t> &list_pt_id,
            std::vector<std::size_t> &list_eqs_id,
            std::vector<std::size_t> &list_eqs_id_without_ghost
            ) const;
    
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
    void mapEqsIDreduced(
            std::size_t mesh_id, const std::vector<std::size_t> &list_pt_id,
            std::vector<std::size_t> &list_eqs_id,
            std::vector<std::size_t> &list_eqs_id_without_ghost
            ) const;

    /**
     * create a local mapping table
     *
     * @param mesh_id       Mesh ID
     * @param pt_id         Vector of point ID
     * @param local_table   Local mapping table
     */
    void createLocalMappingTable(
        std::size_t mesh_id, const std::vector<std::size_t> &list_pt_id,
        DofEquationIdTable &local_table
        ) const;

    /**
     * return DoF from equation ID
     * @param eqs_id    Equation ID
     * @param var_id    Variable ID
     * @param mesh_id   Mesh ID
     * @param pt_id     Point ID
     */
    void mapDoF(std::size_t eqs_id, std::size_t &var_id, std::size_t &mesh_id, std::size_t &pt_id) const;

    /**
     *
     * @param var_id
     * @param mesh_id
     * @return
     */
    const IEquationIdStorage* getPointEquationIdTable(std::size_t var_id, std::size_t mesh_id) const;

    /**
     *
     * @param var_id
     * @param mesh_id
     * @return
     */
    IEquationIdStorage* getPointEquationIdTable(std::size_t var_id, std::size_t mesh_id) { return _map_var2dof[var_id][mesh_id]; }



protected:
    void setTotalDoFs(std::size_t n) {_total_dofs = n;};
    size_t addVariable(std::size_t mesh_id, IEquationIdStorage* address);

private:
    DISALLOW_COPY_AND_ASSIGN(DofEquationIdTable);


private:
    std::vector<std::map<std::size_t, IEquationIdStorage*> > _map_var2dof; // var -> (mesh, address)
    std::map<std::size_t, std::vector<std::size_t> > _map_msh2var; // mesh -> (var)
    std::map<std::size_t, std::vector<std::size_t> > _ghost_pt; // mesh -> ghost points

    std::size_t _total_dofs;
    std::size_t _total_dofs_without_ghosts;
    bool _is_seq_address;

    DofNumberingType::type _numbering_type;
    DofNumberingType::type _local_numbering_type;
};

} //end

#endif //DOFEQUATIONIDTABLE_H_
