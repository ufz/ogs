/*!
   \file InfoMPI.h
   \brief Declaration of class InfoMPI, which provides basic global running time data in MPI

   \author Wenqing Wang
   \version
   \date Nov Dec 2013

   \copyright
   Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef INFOMPI_H_
#define INFOMPI_H_
namespace BaseLib
{
/*!
   \class InfoMPI

   \brief Provide basic MPI data
*/
class InfoMPI
{
    public:
        InfoMPI() {}

        /*!
            set the number of processors and rank

            \param size number of processors
            \param rank ID of current processor
         */
        static inline void setSizeRank(const int size, const int rank)
        {
            _size = size;
            _rank = rank;
        }

        /// get number of processors in use
        static inline int getSize()
        {
            return _size;
        }

        /// get rank, ID of current processors
        static inline int getRank()
        {
            return _rank;
        }

    private:
        /// Rank size, number of processors
        static  int _size;
        /// Rank, ID of current prosessor
        static  int _rank;
};

} // end namespace
#endif

