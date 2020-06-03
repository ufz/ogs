/**
 * \author Norihiro Watanabe
 * \date   2013-04-18
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/Node.h"

template<typename IndexType>struct SteadyDiffusion2DExample1
{
    using LocalMatrixType =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using LocalVectorType = Eigen::VectorXd;

    class LocalAssemblerData
    {
    public:
        void init(MeshLib::Element const& /*unused*/,
                  std::size_t const /*local_matrix_size*/,
                  LocalMatrixType const& localA,
                  LocalVectorType const& localRhs)
        {
            localA_ = &localA;
            localRhs_ = &localRhs;
        }

        void assemble(std::size_t const id,
                      NumLib::LocalToGlobalIndexMap const& dof_table,
                      double const /*t*/, GlobalVector const& /*x*/,
                      GlobalMatrix& /*M*/, GlobalMatrix& K,
                      GlobalVector& b)
        {
            // The local contributions are computed here, usually, but for this
            // particular test all contributions are equal for all elements and are
            // already stored in the localA_ matrix.

            auto const indices = NumLib::getIndices(id, dof_table);
            auto const r_c_indices =
                NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices,
                                                                indices);

            K.add(r_c_indices, *localA_);
            b.add(indices, *localRhs_);
        }

        LocalMatrixType const& getLocalMatrix() const
        {
            return *localA_;
        }

        LocalVectorType const& getLocalVector() const
        {
            return *localRhs_;
        }

    private:
        LocalMatrixType const* localA_ = nullptr;
        LocalVectorType const* localRhs_ = nullptr;
    };

    static
    void initializeLocalData(const MeshLib::Element& e,
            LocalAssemblerData*& data_ptr,
            std::size_t const local_matrix_size,
            SteadyDiffusion2DExample1 const& example)
    {
        data_ptr = new LocalAssemblerData;
        data_ptr->init(e, local_matrix_size, example.localA_, example.localRhs_);
    }

    SteadyDiffusion2DExample1()
        : localA_(4, 4), localRhs_(4)
    {
        msh = MeshLib::MeshGenerator::generateRegularQuadMesh(2.0, mesh_subdivs);
        for (auto* node : msh->getNodes())
        {
            vec_nodeIDs.push_back(node->getID());
        }
        vec_DirichletBC_id.resize(2 * mesh_stride);
        for (std::size_t i = 0; i < mesh_stride; i++)
        {
            // left side
            vec_DirichletBC_id[i] = i * mesh_stride;

            // right side
            vec_DirichletBC_id[mesh_stride + i] = i * mesh_stride + mesh_subdivs;
        }

        // Local assembler matrix and vector are equal for all quad elements.
        {
            localA_(0,0) = 4.0; localA_(0,1) = -1.0; localA_(0,2) = -2.0; localA_(0,3) = -1.0;
            localA_(1,1) = 4.0; localA_(1,2) = -1.0; localA_(1,3) = -2.0;
            localA_(2,2) = 4.0; localA_(2,3) = -1.0;
            localA_(3,3) = 4.0;

            // copy upper triangle to lower to localA for symmetry
            for (std::size_t i = 0; i < 4; i++)
            {
                for (std::size_t j = 0; j < i; j++)
                {
                    localA_(i, j) = localA_(j, i);
                }
            }

            //localA_ *= 1.e-11/6.0;
            for (std::size_t i = 0; i < 4; i++)
            {
                for (std::size_t j = 0; j < 4; j++)
                {
                    localA_(i, j) *= 1.e-11 / 6.0;
                }
            }

            // Fill rhs with zero;
            for (std::size_t i = 0; i < 4; i++)
            {
                localRhs_[i] = 0;
            }
        }

        vec_DirichletBC_value.resize(2 * mesh_stride);
        std::fill_n(vec_DirichletBC_value.begin(), mesh_stride, 0);
        std::fill_n(vec_DirichletBC_value.begin() + mesh_stride, mesh_stride, 1);
        exact_solutions.resize(dim_eqs);
        for (std::size_t i = 0; i < mesh_stride; i++)
        {
            for (std::size_t j = 0; j < mesh_stride; j++)
            {
                exact_solutions[i * mesh_stride + j] = j * 1. / mesh_subdivs;
            }
        }
    }

    ~SteadyDiffusion2DExample1()
    {
        delete msh;
    }

    static const std::size_t mesh_subdivs = 10;
    static const std::size_t mesh_stride = mesh_subdivs + 1;
    static const std::size_t dim_eqs = mesh_stride * mesh_stride;
    MeshLib::Mesh* msh;
    std::vector<IndexType> vec_DirichletBC_id;
    std::vector<double> vec_DirichletBC_value;
    std::vector<double> exact_solutions;
    std::vector<std::size_t> vec_nodeIDs;

    LocalMatrixType localA_;
    LocalVectorType localRhs_;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
