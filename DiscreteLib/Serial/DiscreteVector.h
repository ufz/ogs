/**
 * \file   DiscreteVector.h
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

#include <vector>
#include <cstddef>

#include "BaseLib/CodingTools.h"
#include "MeshLib/Mesh.h"
#include "DiscreteLib/Core/IDiscreteVector.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/Core/IDiscreteVectorAssembler.h"

namespace DiscreteLib
{
class DofEquationIdTable;

/**
 * \brief Vector container for single memory
 */
template<typename T>
class DiscreteVector : public IDiscreteVector<T>
{
protected:
    explicit DiscreteVector(MeshLib::Mesh* msh) : _data(0), _msh(msh), _sys(0) {};
    explicit DiscreteVector(std::size_t n, IDiscreteSystem* sys) : _data(n), _msh(&sys->getMesh()), _sys(sys) {};
    DiscreteVector(MeshLib::Mesh* msh, std::size_t n) : _data(n), _msh(msh), _sys(0) {};

public:
    static DiscreteVector<T>* createInstance(IDiscreteSystem &sys, std::size_t n)
    {
        DiscreteVector<T>* vec = new DiscreteVector<T>(n, &sys);
        sys.addVector(vec);
        return vec;
    }
    virtual ~DiscreteVector() {};

    virtual DiscreteVector<T>* clone() const
    {
        DiscreteVector<T>* vec = DiscreteVector<T>::createInstance(*_sys, _data.size());
        *vec = (*this);
        return vec;
    }

    void resize(std::size_t n) {_data.resize(n);};
    virtual std::size_t size() const {return _data.size();};

    virtual DiscreteVector<T>& operator= (T v)
    {
        for (std::size_t i=getRangeBegin(); i<getRangeEnd(); i++)
            (*this)[i] = v;
        return *this;
    }

    virtual T& operator[] (std::size_t idx)
    {
        return _data[idx];
    }
    virtual const T& operator[] (std::size_t idx) const
    {
        return _data[idx];
    }

    typename std::vector<T>::iterator begin() 
    {
        return _data.begin();
    }
    typename std::vector<T>::iterator end() 
    {
        return _data.end();
    }
    virtual std::size_t getRangeBegin() const
    {
        return 0;
    }
    virtual std::size_t getRangeEnd() const
    {
        return _data.size();
    }

    /// construct
    virtual void construct(const DofEquationIdTable &dofEquationIdTable, const IDiscreteVectorAssembler<T>& assembler)
    {
        assembler.assembly(*_msh, dofEquationIdTable, *this);
    }

protected:
    std::vector<T> _data;
    const MeshLib::Mesh* _msh;
    IDiscreteSystem* _sys;
};

} // end
