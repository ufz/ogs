/*!
   \file InforMPI.h
   \brief Declaration of class InforMPI, which provides basic running time data in MPI

   \author Wenqing Wang
   \version
   \date Nov Dec 2013

   \copyright
   Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef INFORMPI_H_
#define INFORMPI_H_
namespace BaseLib
{
/*!
   \class InforMPI

   \brief Provide basic MPI data
*/
class InforMPI
{
   public:
      InforMPI() {}

      /// get number of computers cores in use
      /*!
          set the number of processors and rank

          \param size number of processors
          \rank ID of current processor
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

