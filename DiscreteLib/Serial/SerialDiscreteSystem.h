/**
 * \file   SerialDiscreteSystem.h
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
#ifndef SERIALDISCRETESYSTEM_H_
#define SERIALDISCRETESYSTEM_H_

#include "BaseLib/CodingTools.h"

#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/DiscreteEnums.h"
#include "SerialDiscreteVector.h"
#include "SerialDiscreteLinearSystem.h"
#include "SequentialElementWiseLinearSystemAssembler.h"
#include "SequentialElementWiseVectorAssembler.h"

namespace MeshLib
{
class Mesh;
}

namespace DiscreteLib
{

class DofEquationIdTable;

/**
 * \brief Non-parallel version of a discrete system
 *
 */
class SerialDiscreteSystem : public IDiscreteSystem
{
public:
    //---------------------------------------------------------------
    // Declaration of sub object types
    //---------------------------------------------------------------
    template<typename T_LINEAR_SOLVER, typename T_SPARSITY_BUILDER>
    struct MyLinearSystem   {
        typedef SerialDiscreteLinearSystem<T_LINEAR_SOLVER, T_SPARSITY_BUILDER> type;
    };

    template<typename T>
    struct MyVector
    {
        typedef SerialDiscreteVector<T> type;
    };

    template <typename T_UPDATER, typename T_LINEAR_SOLVER>
    struct MyLinearSystemAssembler
    {
        typedef SequentialElementWiseLinearSystemAssembler<T_UPDATER, T_LINEAR_SOLVER> type;
    };

    template <typename T_VALUE, typename T_UPDATER>
    struct MyVectorAssembler
    {
        typedef SequentialElementWiseVectorAssembler<T_VALUE, T_UPDATER> type;
    };

    //---------------------------------------------------------------
    //
    //---------------------------------------------------------------
    /**
     * return the system type
     * @return
     */
    static DiscreteSystemType::type getSystemType() { return DiscreteSystemType::Serial;};

    /**
     * Constructor
     * @param msh  pointer to a mesh which represents the discrete system by nodes or elements or edges
     */
    explicit SerialDiscreteSystem(const MeshLib::Mesh* msh) : _msh(msh) {};

    /**
     *
     */
    virtual ~SerialDiscreteSystem()
    {
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
    typename MyLinearSystem<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>::type* createLinearSystem(DofEquationIdTable* dofManager)
    {
        return MyLinearSystem<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>::type::createInstance(*this, _msh, dofManager);
    }

    /// create a new vector
    /// @param n    vector length
    /// @return vector object
    template<typename T>
    typename MyVector<T>::type* createVector(size_t n)
    {
        //return MyVector<T>::type::createInstance(*this, n);
        typedef typename MyVector<T>::type MyVectorType;
        MyVectorType* vec = new MyVectorType(n, this);
        addVector(vec);
        return vec;
    };

private:
    DISALLOW_COPY_AND_ASSIGN(SerialDiscreteSystem);

protected:
    const MeshLib::Mesh* _msh;
};

} //end

#endif //SERIALDISCRETESYSTEM_H_
