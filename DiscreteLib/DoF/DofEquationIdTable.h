/**
 * \file   DofEquationIdTable.h
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief  Helper macros.
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
#include "SequentialEquationIdStorage.h"
#include "RandomEquationIdStorage.h"


namespace DiscreteLib
{

/**
 * \brief DoF and equation index table
 *
 * This class contains a table which relates DoFs and index in a equation.
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
    /// add DoFs
    /// @param var_id       variable type
    /// @param mesh_id      mesh id
    /// @param list_dof_pt_id   a list of discrete point id defined as DoFs
    std::size_t addVariableDoFs(std::size_t mesh_id, const std::vector<std::size_t> &list_dof_pt_id);

    /// add DoFs
    /// @param var_id           variable type
    /// @param mesh_id          mesh id
    /// @param dof_pt_id_begin  the beginning id of discrete points
    /// @param dof_pt_count     the number of discrete points
    std::size_t addVariableDoFs(std::size_t mesh_id, std::size_t dof_pt_id_begin, std::size_t dof_pt_count);

    /// add DoFs
    void addVariableDoFs(std::size_t n_var, std::size_t mesh_id, std::size_t dof_pt_id_begin, std::size_t dof_pt_count);

    /// deactivate DoFs
    void deactivateDoFs(std::size_t var_id, std::size_t mesh_id, std::size_t pt_id);

    /// set ghost points for overlapped regions appeared with domain-decomposition methods
    void setGhostPoints(std::size_t mesh_id, std::vector<std::size_t> &list_pt_id);

    //---------------------------------------------------------------------------------------------
    // construct the table
    //---------------------------------------------------------------------------------------------
    /// numbering DoFs
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
    /// get the number of variables associated with the given mesh
    std::size_t getNumberOfVariables(std::size_t mesh_id) {return _map_msh2var[mesh_id].size();};

    //# Ghost 
    /// is this point a ghost?
    bool isGhostPoint(std::size_t mesh_id, std::size_t pt_id) const;
    /// get the number of ghost points
    std::size_t getNumberOfGhostPoints(std::size_t mesh_id) {return _ghost_pt[mesh_id].size();}
    /// get the total number of ghost points
    std::size_t getTotalNumberOfGhostPoints() const;

    //# DoFs
    /// is this DoF active?
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
     * @param var_id
     * @param mesh_id
     * @param pt_id
     * @return equation id or npos if no match
     */
    std::size_t mapEqsID(std::size_t var_id, std::size_t mesh_id, std::size_t pt_id) const;

    /**
     * return a list of equation id
     *
     * @param var_id
     * @param mesh_id
     * @param pt_id
     * @return equation id or npos if no match
     */
    void mapEqsID(std::size_t var_id, std::size_t mesh_id, const std::vector<std::size_t> &pt_id, std::vector<std::size_t> &eqs_id) const;

    /**
     * return a list of equation id for all variables
     *
     * @param mesh_id
     * @param pt_id
     * @param eqs_id
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
     * @param mesh_id
     * @param pt_id
     * @param eqs_id
     * @param eqs_id_without_ghost
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
     *
     * @param mesh_id
     * @param pt_id
     * @param local_table
     */
    void createLocalMappingTable(
        std::size_t mesh_id, const std::vector<std::size_t> &list_pt_id,
        DofEquationIdTable &local_table
        ) const;

    /**
     *
     * @param eqs_id
     * @param var_id
     * @param mesh_id
     * @param pt_id
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

    /**
     *
     * @param numbering
     */
    void setNumberingType(DofNumberingType::type numbering) {_numbering_type = numbering;};

    /**
     *
     * @return
     */
    DofNumberingType::type getNumberingType() const {return _numbering_type;};

    /**
     *
     * @param numbering
     */
    void setLocalNumberingType(DofNumberingType::type numbering) {_local_numbering_type = numbering;};

    /**
     *
     * @return
     */
    DofNumberingType::type getLocalNumberingType() const {return _local_numbering_type;};

protected:
    void setTotalDoFs(std::size_t n) {_total_dofs = n;};
    //std::map<std::size_t, std::vector<DofMap*> >& getMapMsh2Dof() {return _map_msh2dof;}
    /// get the DoF mapping for the given variable

private:
    DISALLOW_COPY_AND_ASSIGN(DofEquationIdTable);


private:
    std::vector<std::map<std::size_t, IEquationIdStorage*> > _map_var2dof; // var -> (mesh, address)
    //std::vector<Base::IMappedAddress*> _map_var2dof;
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
