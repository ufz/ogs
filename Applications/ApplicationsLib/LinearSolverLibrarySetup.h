// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

#include "BaseLib/ExportSymbol.h"

namespace ApplicationsLib
{
/// The LinearSolverLibrarySetup takes care of proper initialization and
/// shutting down of an external linear solver library. The concrete
/// implementation is chosen by the build system.
/// An object of this class must be created at the beginning of the scope where
/// it is used. When the scope closes (or the object is destroyed explicitly)
/// library shutting down functions are automatically called.
/// The default implementation is empty providing polymorphic behaviour when
/// using this class.
struct LinearSolverLibrarySetup
{
    OGS_EXPORT_SYMBOL std::shared_ptr<LinearSolverLibrarySetup> static create(
        int argc, char* argv[]);

protected:
    LinearSolverLibrarySetup() = default;
    virtual ~LinearSolverLibrarySetup();

private:
    LinearSolverLibrarySetup(const LinearSolverLibrarySetup&) = delete;
    LinearSolverLibrarySetup& operator=(const LinearSolverLibrarySetup&) =
        delete;
};
}  // namespace ApplicationsLib
