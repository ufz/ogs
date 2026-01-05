// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <QObject>
#include <vtkDataSetMapper.h>

/// @brief Simply wraps vtkDataSetMapper as a Qt object to enable slot connections.
class QVtkDataSetMapper : public QObject, public vtkDataSetMapper
{
    Q_OBJECT

public:
    /// @brief Create new objects with New() because of VTKs reference counting.
    static QVtkDataSetMapper* New();

    vtkTypeMacro(QVtkDataSetMapper, vtkDataSetMapper);

    /// @brief Prints information about itself.
    void PrintSelf(ostream& os, vtkIndent indent) override;

public slots:
    /// @brief Sets the scalar visibility on this mapper.
    virtual void SetScalarVisibility(bool on);
    void SetScalarVisibility(int on) override
    {
        SetScalarVisibility(static_cast<bool>(on));
    }

protected:
    /// @brief Constructor.
    QVtkDataSetMapper();

    /// @brief Destructor.
    ~QVtkDataSetMapper() override;

private:
    QVtkDataSetMapper(const QVtkDataSetMapper&); // Not implemented.
    void operator=(const QVtkDataSetMapper&); // Not implemented
};
