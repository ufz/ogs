/**
 * \date   2014-06-11
 * \brief  Implementation of tests.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <fstream>
#include <sstream>

#include <tclap/CmdLine.h>
#include <logog/include/logog.hpp>
#include <logog/include/formatter.hpp>

#include "BaseLib/LogogSimpleFormatter.h"
#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"

int main(int argc, char *argv[])
{
    LOGOG_INITIALIZE();
    BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
    logog::Cout *logogCout(new logog::Cout);
    logogCout->SetFormatter(*custom_format);

    TCLAP::CmdLine cmd("Simple direct matrix solver test.\n\
        It consists of the following steps:\n\
        (1) Read a matrix A from ascii format\n\
        (2) Set all entries of a vector x to one and compute b = A * x\n\
        (3) Solve the system of linear equations -> result have to be (1,...,1)", ' ', "0.1");
    TCLAP::ValueArg<std::string> matrix_arg("m", "matrix", "input matrix file (ascii format)", true, "", "string");
    cmd.add( matrix_arg );
    cmd.parse( argc, argv );

    // *** reading dense matrix in ascii format from file
    std::string const fname_mat(matrix_arg.getValue());
    std::ifstream in(fname_mat.c_str());
    if (!in) {
        INFO("error reading matrix from %s", fname_mat.c_str());
        return -1;
    }
    INFO("reading matrix from %s ...", fname_mat.c_str());
    std::size_t n_rows(0), n_cols(0);
    in >> n_rows;
    in >> n_cols;
    MathLib::DenseMatrix<double, std::size_t> mat(n_rows, n_cols);
    for (std::size_t i(0); i<mat.getNumberOfRows(); ++i) {
        for (std::size_t j(0); j<mat.getNumberOfColumns(); ++j) {
            in >> mat(i,j);
        }
    }
    {
        std::stringstream stream;
        stream << mat;
        INFO("read matrix:\n%s", stream.str().c_str());
    }

    std::vector<double> x(n_cols,1.0), b;
    b.resize(n_rows);
    b = mat * x;

    MathLib::GaussAlgorithm<MathLib::DenseMatrix<double, std::size_t>> gauss;
    gauss.solve(mat, b, true);

    {
        std::stringstream stream;
        std::copy(b.begin(), b.end(), std::ostream_iterator<double>(stream, " "));
        stream << std::endl;
        INFO("solution vector:\n%s", stream.str().c_str());
    }

    delete custom_format;
    delete logogCout;
    LOGOG_SHUTDOWN();

    return 0;
}

