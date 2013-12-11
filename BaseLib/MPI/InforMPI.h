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
      static inline void setSizeRank(const int size, const int rank)
      {
         _size = size;
         _rank = rank;
      }

      /// get number of computers cores in use
      static inline int getSize()
      {
         return _size;
      }

      /// get rank
      static inline int getRank()
      {
         return _rank;
      }

   private:

      /// Rank size
      static  int _size;
      /// Rank
      static  int _rank;

};

} // end namespace
#endif

