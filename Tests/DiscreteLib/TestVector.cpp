/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief  Implementation tests.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "BaseLib/CodingTools.h"
#include "MeshLib/MeshGenerator.h"
#include "DiscreteLib/ElementWiseManipulator/IElemenetWiseVectorLocalAssembler.h"
#include "DiscreteLib/ElementWiseManipulator/ElementWiseVectorUpdater.h"
#include "DiscreteLib/Serial/SerialDiscreteSystem.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace MathLib;
using namespace MeshLib;
using namespace DiscreteLib;

namespace
{
inline void ASSERT_DOUBLE_ARRAY_EQ(const double* Expected, const double* Actual, size_t N, double epsilon=1.0e-8) {
    for (size_t i=0; i<N; i++) \
        ASSERT_NEAR(Expected[i], Actual[i], epsilon);
}

class LocalVectorAssembler1 : public DiscreteLib::IElemenetWiseVectorLocalAssembler
{
public:
    virtual ~LocalVectorAssembler1() {};
    virtual void assembly(const MeshLib::Element &/*e*/, LocalVector &local_v)
    {
        for (std::size_t i=0; i<4; i++)
            local_v[i] = i;
    }
};

} //namespace


TEST(DiscreteLib, VectorSerial1)
{
    std::unique_ptr<MeshLib::Mesh> msh(MeshGenerator::generateRegularQuadMesh(2.0, 2, .0, .0, .0));

    SerialDiscreteSystem dis(msh.get());
    SerialDiscreteVector<double> *v = dis.createVector<double>(msh->getNNodes());

    double expected[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    for (std::size_t i=0; i<9; i++)
        (*v)[i] = expected[i];
    
    ASSERT_EQ(9u, v->size());
    ASSERT_DOUBLE_ARRAY_EQ(expected, &(*v)[0], 9);
}


TEST(DiscreteLib, VectorSerial2)
{
    std::unique_ptr<MeshLib::Mesh> msh(MeshGenerator::generateRegularQuadMesh(2.0, 2, .0, .0, .0));

    SerialDiscreteSystem dis(msh.get());
    SerialDiscreteVector<double> *v = dis.createVector<double>(msh->getNNodes());


    // set up a DoF table
    DofEquationIdTable dofMappingTable;
    dofMappingTable.addVariableDoFs(msh->getID(), 0, msh->getNNodes());
    dofMappingTable.construct(DofNumberingType::BY_VARIABLE);

    LocalVectorAssembler1 local_assembler;
    typedef DiscreteLib::ElementWiseVectorUpdater<double, LocalVectorAssembler1> MyElementWiseVectorUpdater;
    MyElementWiseVectorUpdater updater(msh->getID(), local_assembler);
    SequentialElementWiseVectorAssembler<double, MyElementWiseVectorUpdater> e_assembler(updater);
    v->construct(dofMappingTable, e_assembler);


    double expected[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    for (std::size_t i=0; i<9; i++)
        (*v)[i] = expected[i];

    ASSERT_EQ(9u, v->size());
    ASSERT_DOUBLE_ARRAY_EQ(expected, &(*v)[0], 9);
}

