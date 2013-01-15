/**
 * \file   DiscreteSystem.h
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
#pragma once

#include "BaseLib/CodingTools.h"

#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/DiscreteEnums.h"
#include "DiscreteVector.h"
#include "DiscreteLinearEquation.h"
#include "SequentialElementWiseLinearEquationAssembler.h"
#include "SequentialElementWiseVectorAssembler.h"

namespace MeshLib
{
class Mesh;
}

namespace DiscreteLib
{

class DofEquationIdTable;

/**
 * \brief Discrete system based on a single mesh
 *
 * - Mesh
 */
class DiscreteSystem : public IDiscreteSystem
{
public:
    template<typename T_LINEAR_SOLVER, typename T_SPARSITY_BUILDER>
    struct MyLinearEquation
    {
        typedef DiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER> type;
    };

    template<typename T>
    struct MyVector
    {
        typedef DiscreteVector<T> type;
    };

    template <typename T_UPDATER, typename T_LINEAR_SOLVER>
    struct MyLinearEquationAssembler
    {
        typedef SequentialElementWiseLinearEquationAssembler<T_UPDATER, T_LINEAR_SOLVER> type;
    };

    template <typename T_VALUE, typename T_UPDATER>
    struct MyVectorAssembler
    {
        typedef SequentialElementWiseVectorAssembler<T_VALUE, T_UPDATER> type;
    };

    /**
     * 
     * @return
     */
    static DiscreteSystemType::type getSystemType() { return DiscreteSystemType::Serial;};

    /// constructor
    /// \param msh  a mesh to represent discrete systems by nodes or elements or edges
    explicit DiscreteSystem(const MeshLib::Mesh* msh) : _msh(msh) {};

    ///
    virtual ~DiscreteSystem()
    {
        //BaseLib::releaseObject(_msh);
    };

    /// return this mesh
    virtual const MeshLib::Mesh& getMesh() const {return *_msh;};

    /// create a new linear equation
    /// @tparam T_LINEAR_EQUATION
    /// @tparam T_LINEAR_SOLVER
    /// @tparam T_SPARSITY_BUILDER
    /// @param linear_solver         Linear solver
    /// @param dofManager            Equation index table
    template <class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
    typename MyLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>::type* createLinearEquation(T_LINEAR_SOLVER* linear_solver, DofEquationIdTable* dofManager)
    {
        return MyLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>::type::createInstance(*this, _msh, linear_solver, dofManager);
    }

    /// create a new vector
    /// @param n    vector length
    /// @return vector object
    template<typename T>
    typename MyVector<T>::type* createVector(size_t n)
    {
        return MyVector<T>::type::createInstance(*this, n);
    };

private:
    DISALLOW_COPY_AND_ASSIGN(DiscreteSystem);

protected:
    const MeshLib::Mesh* _msh;
};

} //end
