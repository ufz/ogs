/**
 * \file   SerialDiscreteVector.h
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
#ifndef SERIALDISCRETEVECTOR_H_
#define SERIALDISCRETEVECTOR_H_

#include <vector>

#include "BaseLib/CodingTools.h"
#include "DiscreteLib/Core/IDiscreteVector.h"
#include "DiscreteLib/Core/IDiscreteVectorAssembler.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"

namespace MeshLib
{
class Mesh;
}

namespace DiscreteLib
{

//forward declaration
class DofEquationIdTable;
class SerialDiscreteSystem;

/**
 * \brief Non-parallelized version of vector class
 */
template<typename T>
class SerialDiscreteVector : public IDiscreteVector<T>
{
public:
    // only SerialDiscreteSystem can call constructor
    friend class SerialDiscreteSystem;

    //---------------------------------------------------------------
    // override member function
    //---------------------------------------------------------------
    /**
     *
     */
    virtual ~SerialDiscreteVector() {};

    /**
     * return a clone of this object
     * @return
     */
    virtual SerialDiscreteVector<T>* clone() const
    {
        SerialDiscreteVector<T>* vec = new SerialDiscreteVector<T>(*this);
        _dis_resource->addVector(vec);
        return vec;
    }

    /**
     * return the vector size
     * @return
     */
    virtual std::size_t size() const {return _data.size();};

    /**
     * operator =
     * @param v
     * @return
     */
    virtual SerialDiscreteVector<T>& operator= (T v)
    {
        for (std::size_t i=getRangeBegin(); i<getRangeEnd(); i++)
            (*this)[i] = v;
        return *this;
    }

    /**
     * access an entry
     * @param idx
     * @return
     */
    virtual T& operator[] (std::size_t idx)
    {
        return _data[idx];
    }

    /**
     * access an entry
     * @param idx
     * @return
     */
    virtual const T& operator[] (std::size_t idx) const
    {
        return _data[idx];
    }

    /**
     * return the starting index of the valid range
     * @return
     */
    virtual std::size_t getRangeBegin() const
    {
        return 0;
    }

    /**
     * return the end index of the valid range
     * @return
     */
    virtual std::size_t getRangeEnd() const
    {
        return _data.size();
    }

    /**
     * construct the vector entries
     *
     * @param dofEquationIdTable    A mapping table between DoFs and equation index
     * @param assembler             Assembler
     */
    virtual void construct(const DofEquationIdTable &dofEquationIdTable, IDiscreteVectorAssembler<T>& assembler)
    {
        assembler.assembly(*_msh, dofEquationIdTable, *this);
    }

protected:
    /**
     * Constructor
     * @param n             Vector size
     * @param msh           Mesh object
     * @param dis_resource  Discrete data resource manager
     */
    SerialDiscreteVector(std::size_t n, const MeshLib::Mesh* msh, DiscreteResourceManager* dis_resource)
    : _data(n), _msh(msh), _dis_resource(dis_resource) {};

    /**
     * Copy constructor
     * @param src
     */
    SerialDiscreteVector(const SerialDiscreteVector &src)
    : _data(src._data), _msh(src._msh), _dis_resource(src._dis_resource) {};

protected:
    std::vector<T> _data;
    const MeshLib::Mesh* _msh;
    DiscreteResourceManager* _dis_resource;
};

} // end

#endif //SERIALDISCRETEVECTOR_H_
