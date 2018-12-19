/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
struct SigmaIntegrationPointWriter final : public IntegrationPointWriter
{
    SigmaIntegrationPointWriter(
        int const n_components,
        int const integration_order,
        std::function<std::vector<std::vector<double>>()>
            callback)
        : _callback(callback),
          _integration_order(integration_order),
          _n_components(n_components)
    {
    }

    int numberOfComponents() const override { return _n_components; }
    int integrationOrder() const override { return _integration_order; }

    std::string name() const override
    {
        // TODO (naumov) remove ip suffix. Probably needs modification of the
        // mesh properties, s.t. there is no "overlapping" with cell/point data.
        // See getOrCreateMeshProperty.
        return "sigma_ip";
    }

    std::vector<std::vector<double>> values() const override
    {
        return _callback();
    }

private:
    std::function<std::vector<std::vector<double>>()> _callback;
    int const _integration_order;
    int const _n_components;
};

struct KappaDIntegrationPointWriter final : public IntegrationPointWriter
{
    KappaDIntegrationPointWriter(
        int const integration_order,
        std::function<std::vector<std::vector<double>>()>
            callback)
        : _callback(callback), _integration_order(integration_order)
    {
    }

    int numberOfComponents() const override { return 1; }
    int integrationOrder() const override { return _integration_order; }

    std::string name() const override
    {
        // TODO (naumov) remove ip suffix. Probably needs modification of the
        // mesh properties, s.t. there is no "overlapping" with cell/point data.
        // See getOrCreateMeshProperty.
        return "kappa_d_ip";
    }

    std::vector<std::vector<double>> values() const override
    {
        return _callback();
    }

private:
    std::function<std::vector<std::vector<double>>()> _callback;
    int const _integration_order;
};

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
