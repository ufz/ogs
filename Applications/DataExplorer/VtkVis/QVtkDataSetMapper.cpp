// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// ** INCLUDES **
#include "QVtkDataSetMapper.h"

#include <vtkObjectFactory.h>

vtkStandardNewMacro(QVtkDataSetMapper);

QVtkDataSetMapper::QVtkDataSetMapper() : QObject(nullptr) {}

QVtkDataSetMapper::~QVtkDataSetMapper() = default;

void QVtkDataSetMapper::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}

void QVtkDataSetMapper::SetScalarVisibility(bool on)
{
    vtkDataSetMapper::SetScalarVisibility(static_cast<int>(on));
}
