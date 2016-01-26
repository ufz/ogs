/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
#define PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_

#include <cassert>

#include <boost/optional.hpp>

#include "AssemblerLib/LocalAssemblerBuilder.h"
#include "AssemblerLib/LocalDataInitializer.h"

#include "FileIO/VtkIO/VtuInterface.h"

#include "GroundwaterFlowFEM.h"
#include "Parameter.h"
#include "Process.h"

namespace MeshLib
{
    class Element;
    class Mesh;
    template <typename PROP_VAL_TYPE> class PropertyVector;
}

namespace ProcessLib
{

template<typename GlobalSetup>
class GroundwaterFlowProcess : public Process<GlobalSetup>
{
public:
    GroundwaterFlowProcess(MeshLib::Mesh& mesh,
            std::vector<ProcessVariable> const& variables,
            std::vector<std::unique_ptr<ParameterBase>> const& parameters,
            BaseLib::ConfigTreeNew const& config)
        : Process<GlobalSetup>(mesh)
    {
        config.checkConfParam("type", "GROUNDWATER_FLOW");

        DBUG("Create GroundwaterFlowProcess.");

        // Process variable.
        {
            // Find the corresponding process variable.
            std::string const name = config.getConfParam<std::string>("process_variable");

            auto variable = std::find_if(variables.cbegin(), variables.cend(),
                    [&name](ProcessVariable const& v) {
                        return v.getName() == name;
                    });

            if (variable == variables.end())
                ERR("Expected process variable \'%s\' not found in provided variables list.",
                    name.c_str());

            DBUG("Associate hydraulic_head with process variable \'%s\'.",
                name.c_str());
            this->_process_variables.emplace_back(
                const_cast<ProcessVariable*>(&*variable));
        }

        // Hydraulic conductivity parameter.
        {
            // find hydraulic_conductivity in process config
            auto const name = config.getConfParam<std::string>("hydraulic_conductivity");

            // find corresponding parameter by name
            auto const parameter =
                std::find_if(parameters.cbegin(), parameters.cend(),
                             [&name](std::unique_ptr<ParameterBase> const& p)
                             {
                                 return p->name == name;
                             });

            if (parameter == parameters.end())
            {
                ERR("Could not find required parameter config for \'%s\' "
                    "among read parameters.",
                    name.c_str());
                std::abort();
            }

            _hydraulic_conductivity =
                dynamic_cast<const Parameter<double, const MeshLib::Element&>*>(
                    parameter->get());
            if (!_hydraulic_conductivity)
            {
                ERR("The hydraulic conductivity parameter is of incompatible "
                    "type.");
                std::abort();
            }
        }

        // Linear solver options
        if (auto linear_solver_options =
                config.getConfSubtreeOptional("linear_solver"))
            Process<GlobalSetup>::setLinearSolverOptions(
                std::move(*linear_solver_options));
    }

    template <unsigned GlobalDim>
    void createLocalAssemblers()
    {
        DBUG("Create local assemblers.");
        // Populate the vector of local assemblers.
        _local_assemblers.resize(this->_mesh.getNElements());
        // Shape matrices initializer
        using LocalDataInitializer = AssemblerLib::LocalDataInitializer<
            GroundwaterFlow::LocalAssemblerDataInterface,
            GroundwaterFlow::LocalAssemblerData,
            typename GlobalSetup::MatrixType,
            typename GlobalSetup::VectorType,
            GlobalDim>;

        LocalDataInitializer initializer;

        using LocalAssemblerBuilder =
            AssemblerLib::LocalAssemblerBuilder<
                MeshLib::Element,
                LocalDataInitializer>;

        LocalAssemblerBuilder local_asm_builder(
            initializer, *this->_local_to_global_index_map);

        DBUG("Calling local assembler builder for all mesh elements.");
        this->_global_setup.execute(
                local_asm_builder,
                this->_mesh.getElements(),
                _local_assemblers,
                *_hydraulic_conductivity,
                this->_integration_order);
    }

    std::string getLinearSolverName() const override
    {
        return "gw_";
    }

    void createLocalAssemblers() override
    {
        if (this->_mesh.getDimension()==1)
            createLocalAssemblers<1>();
        else if (this->_mesh.getDimension()==2)
            createLocalAssemblers<2>();
        else if (this->_mesh.getDimension()==3)
            createLocalAssemblers<3>();
        else
            assert(false);
    }

    bool assemble(const double /*delta_t*/) override
    {
        DBUG("Assemble GroundwaterFlowProcess.");

        *this->_rhs = 0;   // This resets the whole vector.

        // Call global assembler for each local assembly item.
        this->_global_setup.execute(*this->_global_assembler,
                                    _local_assemblers);

        return true;
    }

    void post(std::string const& file_name) override
    {
        DBUG("Postprocessing GroundwaterFlowProcess.");

        std::string const property_name = "Result";

        // Get or create a property vector for results.
        boost::optional<MeshLib::PropertyVector<double>&> result;
        if (this->_mesh.getProperties().hasPropertyVector(property_name))
        {
            result = this->_mesh.getProperties().template
                getPropertyVector<double>(property_name);
        }
        else
        {
            result = this->_mesh.getProperties().template
                createNewPropertyVector<double>(property_name,
                    MeshLib::MeshItemType::Node);
            result->resize(this->_x->size());
        }
        assert(result && result->size() == this->_x->size());

        // Copy result
#ifdef USE_PETSC
        double* loc_x = _x->getLocalVector();
        for (std::size_t i=0; i<_x->getLocalSize() + _x->getGhostSize(); i++)
        {
            (*result)[i] = loc_x[i];
        }
        _x->restoreArray(loc_x);
#else
        for (std::size_t i = 0; i < this->_x->size(); ++i)
            (*result)[i] = (*this->_x)[i];
#endif

        // Write output file
        FileIO::VtuInterface vtu_interface(&this->_mesh, vtkXMLWriter::Binary, true);
        vtu_interface.writeToFile(file_name);
    }

    void postTimestep(std::string const& file_name, const unsigned /*timestep*/) override
    {
        post(file_name);
    }

    ~GroundwaterFlowProcess()
    {
        for (auto p : _local_assemblers)
            delete p;
    }

private:
    Parameter<double, MeshLib::Element const&> const* _hydraulic_conductivity = nullptr;

    using LocalAssembler = GroundwaterFlow::LocalAssemblerDataInterface<
        typename GlobalSetup::MatrixType, typename GlobalSetup::VectorType>;

    std::vector<LocalAssembler*> _local_assemblers;
};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
