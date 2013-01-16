/**
 * \file   IDiscreteVector.h
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

#ifndef IDISCRETEVECTOR_H_
#define IDISCRETEVECTOR_H_

#include <vector>

#include "IDiscreteVectorBase.h"

namespace DiscreteLib
{

// forward declaration
template <typename T> class IDiscreteVectorAssembler;
class DofEquationIdTable;

/**
 * \brief Interface of discrete data with specific data type
 * 
 * \tparam T    Data type, e.g. int, double
 */
template <typename T>
class IDiscreteVector : public IDiscreteVectorBase
{
public:
    typedef IDiscreteVector<T> MyVectorType;

    /// Index indicating out of the range
    static const std::size_t index_npos = -1;
    
    /**
     *
     */
    virtual ~IDiscreteVector() {};

    /**
     * clone this object
     * @return a pointer to the cloned object
     */
    virtual MyVectorType* clone() const = 0;

    /**
     *
     * @return the number of entries
     */
    virtual std::size_t size() const = 0;

    /**
     *
     * @return a start index of the active data range
     */
    virtual std::size_t getRangeBegin() const = 0;

    /**
     *
     * @return an end index of the active data range
     */
    virtual std::size_t getRangeEnd() const = 0;

    /**
     * access data
     * @param idx   index
     * @return value
     */
    virtual T& operator[] (std::size_t idx) = 0;

    /**
     * access data
     * @param idx   index
     * @return value
     */
    virtual const T& operator[] (std::size_t idx) const = 0;

    /**
     * vector operation: set data
     * @param src
     * @return
     */
    virtual MyVectorType& operator= (const MyVectorType &src);

    /**
     * vector operation: add
     * @param v
     */
    virtual void operator+= (const MyVectorType& v);

    /**
     * vector operation: subtract
     * @param v
     */
    virtual void operator-= (const MyVectorType& v);

    /**
     * set all values in this vector
     * @param v
     * @return
     */
    virtual MyVectorType& operator= (T v);

    /**
     * add values to given entries
     * @param vec_pos
     * @param local_v
     */
    virtual void addSubvector(const std::vector<std::size_t> &vec_pos, T* local_v);

    /**
     * set values to given entries
     * @param vec_pos
     * @param v
     */
    virtual void setSubvector(const std::vector<std::size_t> &vec_pos, T v);

    /**
     * construct the vector
     *
     * @param dofEquationIdTable    DoF mapping table
     * @param assembler             Assembler
     */
    virtual void construct( const DofEquationIdTable &dofEquationIdTable,
                            IDiscreteVectorAssembler<T>& assembler) = 0;

};

} // end

//include implementation of template functions
#include "IDiscreteVector.tpp"

#endif //IDISCRETEVECTOR_H_
