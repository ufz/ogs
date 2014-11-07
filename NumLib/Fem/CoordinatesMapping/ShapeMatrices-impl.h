/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


namespace NumLib
{

namespace detail
{
/*
 * zero reset functions for Eigen
 *
 */
template<class T>
void setMatrixZero(T &mat)
{
    mat.setZero(mat.rows(), mat.cols());
}

template<class T>
void setVectorZero(T &vec)
{
    vec.setZero(vec.size());
}

/*
 * "tag dispatching" technique is used below to emulate explicit specialization
 * of a template function in a template class
 */
template <ShapeMatrixType FIELD_TYPE> struct ShapeDataFieldType {};

template <class T_N, class T_DN, class T_J>
inline void setZero(ShapeMatrices<T_N, T_DN, T_J> &shape, ShapeDataFieldType<ShapeMatrixType::N>)
{
    setVectorZero(shape.N);
}

template <class T_N, class T_DN, class T_J>
inline void setZero(ShapeMatrices<T_N, T_DN, T_J> &shape, ShapeDataFieldType<ShapeMatrixType::DNDR>)
{
    setMatrixZero(shape.dNdr);
}

template <class T_N, class T_DN, class T_J>
inline void setZero(ShapeMatrices<T_N, T_DN, T_J> &shape, ShapeDataFieldType<ShapeMatrixType::DNDR_J>)
{
    setZero(shape, ShapeDataFieldType<ShapeMatrixType::DNDR>());
    setMatrixZero(shape.J);
    shape.detJ = .0;
}

template <class T_N, class T_DN, class T_J>
inline void setZero(ShapeMatrices<T_N, T_DN, T_J> &shape, ShapeDataFieldType<ShapeMatrixType::N_J>)
{
    setZero(shape, ShapeDataFieldType<ShapeMatrixType::N>());
    setZero(shape, ShapeDataFieldType<ShapeMatrixType::DNDR_J>());
}

template <class T_N, class T_DN, class T_J>
inline void setZero(ShapeMatrices<T_N, T_DN, T_J> &shape, ShapeDataFieldType<ShapeMatrixType::DNDX>)
{
    setZero(shape, ShapeDataFieldType<ShapeMatrixType::DNDR_J>());
    setMatrixZero(shape.invJ);
    setMatrixZero(shape.dNdx);
}

template <class T_N, class T_DN, class T_J>
inline void setZero(ShapeMatrices<T_N, T_DN, T_J> &shape, ShapeDataFieldType<ShapeMatrixType::ALL>)
{
    setZero(shape, ShapeDataFieldType<ShapeMatrixType::N>());
    setZero(shape, ShapeDataFieldType<ShapeMatrixType::DNDX>());
}

} //detail

template <class T_N, class T_DN, class T_J>
inline void ShapeMatrices<T_N, T_DN, T_J>::setZero()
{
    detail::setZero(*this, detail::ShapeDataFieldType<ShapeMatrixType::ALL>());
}

template <class T_N, class T_DN, class T_J>
template <ShapeMatrixType T_SHAPE_MATRIX_TYPE>
inline void ShapeMatrices<T_N, T_DN, T_J>::setZero()
{
    detail::setZero(*this, detail::ShapeDataFieldType<T_SHAPE_MATRIX_TYPE>());
}

template <class T_N, class T_DN, class T_J>
void ShapeMatrices<T_N, T_DN, T_J>::resize(std::size_t const dim, std::size_t n_nodes)
{
    N.resize(n_nodes);
    dNdr.resize(dim, n_nodes);
    J.resize(dim, dim);
    invJ.resize(dim, dim);
    dNdx.resize(dim, n_nodes);

    setZero();
}

template <class T_N, class T_DN, class T_J>
void ShapeMatrices<T_N, T_DN, T_J>::write(std::ostream& out) const
{
    out << "N   :\n" << N << "\n";
    out << "dNdr:\n" << dNdr << "\n";
    out << "J   :\n" << J << "\n";
    out << "|J| : " << detJ << "\n";
    out << "invJ:\n" << invJ << "\n";
    out << "dNdx:\n" << dNdx << "\n";
}

template <class T_N, class T_DN, class T_J>
std::ostream& operator<< (std::ostream &os, const ShapeMatrices<T_N, T_DN, T_J> &shape)
{
    shape.write (os);
    return os;
}

} // NumLib

