/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EIGENVECTOR_H_
#define EIGENVECTOR_H_

#include <string>
#include <vector>
#include <fstream>

#include <Eigen/Eigen>

namespace MathLib
{

/**
 * Global vector based on Eigen vector
 */
class EigenVector
{
public:
    typedef Eigen::VectorXd RawVectorType;

    /**
     * Constructor for initialization of the number of rows
     * @param length number of rows
     */
    explicit EigenVector(std::size_t length) : _vec(length) {}

    /// copy constructor
    EigenVector(EigenVector const &src) : _vec(src._vec) {}

    /**
     *
     */
    virtual ~EigenVector() {}

    /// return a vector length
    std::size_t size() const { return _vec.size(); }

    /// return a start index of the active data range
    std::size_t getRangeBegin() const { return 0;}

    /// return an end index of the active data range
    std::size_t getRangeEnd() const { return this->size(); }

    /// set all values in this vector
    EigenVector& operator= (double v) { _vec.setConstant(v); return *this; }

    /// set all values in this vector
    EigenVector& operator*= (double v) { _vec *= v; return *this; }

    /// access entry
    double operator[] (std::size_t rowId) const { return get(rowId); }

    /// get entry
    double get(std::size_t rowId) const
    {
        return _vec[rowId];
    }

    /// set entry
    void set(std::size_t rowId, double v)
    {
        _vec[rowId] = v;
    }

    /// add entry
    void add(std::size_t rowId, double v)
    {
        _vec[rowId] += v;
    }

    /// add entries
    template<class T_SUBVEC>
    void add(const std::vector<std::size_t> &pos, const T_SUBVEC &sub_vec)
    {
        for (std::size_t i=0; i<pos.size(); ++i) {
            this->add(pos[i], sub_vec[i]);
        }
    }

    /// printout this equation for debugging
    void write (const std::string &filename) const { std::ofstream os(filename); os << _vec; }

    /// return a raw Lis vector object
    RawVectorType& getRawVector() {return _vec; }

    /// return a raw Lis vector object
    const RawVectorType& getRawVector() const {return _vec; }

    /// vector operation: set data
    EigenVector& operator= (const EigenVector &src) { _vec = static_cast<const EigenVector&>(src)._vec; return *this; }

    /// vector operation: add
    void operator+= (const EigenVector& v) { _vec += static_cast<const EigenVector&>(v)._vec; }

    /// vector operation: subtract
    void operator-= (const EigenVector& v) { _vec -= static_cast<const EigenVector&>(v)._vec; }

private:
    RawVectorType _vec;
};

} // MathLib

#endif //EIGENVECTOR_H_

