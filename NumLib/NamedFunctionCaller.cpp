/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NamedFunctionCaller.h"

#include <algorithm>

#include "BaseLib/uniqueInsert.h"
#include "BaseLib/Error.h"

/// To understand the working of the algorithm the following two lemmata are
/// useful:
/// If a directed graph has a topological order, then it is a Directed Acyclic
/// Graph (DAG).
/// If a graph is a DAG, then it has a node with no incoming edges.
///
/// \note Negative indices in the graph adjacency list are leaves of the graph.
///
/// \param graph represents directed graph given as an adjacency list.
/// \return true if the given directed graph is acyclic.
bool hasTopologicalOrdering(std::vector<std::vector<int>> const& graph)
{
    // Number of incoming edges for each node of the graph
    std::vector<unsigned> number_incoming_edges(graph.size());

    // Count all incoming edges (i,j) for each node (i).
    for (auto const& node_i_adjacencies : graph)
    {
        for (int node_j : node_i_adjacencies)
        {
            if (node_j >= 0)    // ignore negative indices.
                ++number_incoming_edges[node_j];
        }
    }

    // Working queue: a set of nodes with no incoming edges.
    std::vector<std::size_t> q;
    for (std::size_t node_i = 0; node_i < number_incoming_edges.size();
         ++node_i)
    {
        if (number_incoming_edges[node_i] == 0)
        {
            q.push_back(node_i);
        }
    }


    auto num_dependent = number_incoming_edges.size() - q.size();

    while (!q.empty())
    {
        auto const node_i = q.back();
        q.pop_back();

        // Decrement counts for all edges (i,j).
        for (int node_j : graph[node_i])
        {
            if (node_j < 0)
                continue;  // ignore negative indices
            if (--number_incoming_edges[node_j] == 0)
            {  // Add a node without incoming edges to the queue.
                q.push_back(node_j);
                --num_dependent;
            }
        }
    }

    return num_dependent == 0;
}

enum class TraversePosition { StartNode, BetweenChildren, EndNode };

/*! Traverses the graph given by the adjacency list \c map_sink_source in a
 * depth-first manner starting at the node \c sink_fct.
 *
 * At the beginning and end of each node as well as between every two child
 * nodes the given \c callback is called.
 */
template <typename Callback>
void traverse(std::vector<std::vector<int>> const& map_sink_source,
              int sink_fct, Callback&& callback)
{
    assert(sink_fct < static_cast<int>(map_sink_source.size()));
    callback(sink_fct, TraversePosition::StartNode);

    if (sink_fct < 0)
        return;

    auto const& si_so = map_sink_source[sink_fct];
    std::size_t const num_args = si_so.size();
    for (std::size_t sink_arg = 0; sink_arg != num_args; ++sink_arg) {
        if (sink_arg != 0)
            callback(sink_fct, TraversePosition::BetweenChildren);
        traverse(map_sink_source, si_so[sink_arg], callback);
    }

    callback(sink_fct, TraversePosition::EndNode);
}

namespace NumLib
{
NamedFunctionCaller::NamedFunctionCaller(
    std::initializer_list<std::string> unbound_argument_names)
    : _uninitialized(-1 - unbound_argument_names.size())
{
    int idx = -1;
    for (auto arg : unbound_argument_names) {
        BaseLib::insertIfKeyUniqueElseError(
            _map_name_idx, arg, idx,
            "The name of the unbound argument is not unique.");
        --idx;
    }
}

void NamedFunctionCaller::addNamedFunction(NamedFunction&& fct)
{
    DBUG("Adding named function `%s'", fct.getName().c_str());

    BaseLib::insertIfKeyUniqueElseError(
        _map_name_idx, fct.getName(), _named_functions.size(),
        "The name of the function is not unique.");

    _map_sink_source.emplace_back(fct.getArgumentNames().size(),
                                  _uninitialized);
    _named_functions.push_back(std::move(fct));
}

void NamedFunctionCaller::plug(const std::string& sink_fct,
                               const std::string& sink_arg,
                               const std::string& source_fct)
{
    _deferred_plugs.push_back({sink_fct, sink_arg, source_fct});
}

void NamedFunctionCaller::applyPlugs()
{
    while (!_deferred_plugs.empty())
    {
        auto const& plug = _deferred_plugs.back();
        auto const& sink_fct = plug.sink_fct;
        auto const& sink_arg = plug.sink_arg;
        auto const& source = plug.source;

        auto const source_it = _map_name_idx.find(source);
        if (source_it == _map_name_idx.end())
        {
            OGS_FATAL("A function with the name `%s' has not been found.",
                      source.c_str());
        }
        auto const source_idx = source_it->second;

        auto const sink_it = _map_name_idx.find(sink_fct);
        if (sink_it == _map_name_idx.end())
        {
            OGS_FATAL("A function with the name `%s' has not been found.",
                      sink_fct.c_str());
        }
        auto const sink_fct_idx = sink_it->second;

        auto const& sink_args =
            _named_functions[sink_it->second].getArgumentNames();
        auto const sink_arg_it =
            std::find(sink_args.begin(), sink_args.end(), sink_arg);
        if (sink_arg_it == sink_args.end())
        {
            OGS_FATAL(
                "An argument with the name `%s' has not been found for the "
                "function `%s'.",
                sink_arg.c_str(), sink_fct.c_str());
        }
        std::size_t const sink_arg_idx =
            std::distance(sink_args.begin(), sink_arg_it);

        auto& sis_sos = _map_sink_source[sink_fct_idx];
        if (sis_sos[sink_arg_idx] != _uninitialized)
        {
            OGS_FATAL("A dependency for `%s'.`%s' has already been introduced.",
                      sink_fct.c_str(), sink_arg.c_str());
        }
        sis_sos[sink_arg_idx] = source_idx;
        if (!hasTopologicalOrdering(_map_sink_source))
        {
            OGS_FATAL(
                "The call graph being plugged together must be an acyclic "
                "graph. The added dependency for `%s'.`%s' introduces a cycle "
                "into the graph.",
                sink_fct.c_str(), sink_arg.c_str());
        }

        _deferred_plugs.pop_back();
    }
}

double NamedFunctionCaller::call(
    std::size_t function_idx,
    const std::vector<double>& unbound_arguments) const
{
    assert(_deferred_plugs.empty() &&
           "You must call applyPlugs() before this method!");

    auto const& sis_sos = _map_sink_source[function_idx];
    assert(sis_sos.size() ==
           _named_functions[function_idx].getArgumentNames().size());
    std::vector<double> fct_args(sis_sos.size());

    for (std::size_t sink=0; sink<sis_sos.size(); ++sink)
    {
        auto const source = sis_sos[sink];

        if (source >= 0) {
            fct_args[sink] = call(source, unbound_arguments);
        } else {
            assert(source != _uninitialized);
            fct_args[sink] = unbound_arguments[-source-1];
        }
    }

    return _named_functions[function_idx].call(fct_args);
}

std::string NamedFunctionCaller::getCallExpression(
    std::string const& function_name) const
{
    auto const fct_it = _map_name_idx.find(function_name);
    if (fct_it == _map_name_idx.end()) {
        OGS_FATAL("A function with the name `%s' has not been found.",
                  function_name.c_str());
    }

    std::string expr;
    auto callback = [&](int fct_idx, TraversePosition pos)
    {
        switch (pos) {
        case TraversePosition::StartNode:
        {
            if (fct_idx < 0) {
                auto it = std::find_if(
                    _map_name_idx.begin(), _map_name_idx.end(),
                    [fct_idx](std::pair<std::string, int> const& e) {
                        return e.second == fct_idx;
                    });
                if (it == _map_name_idx.end()) {
                    OGS_FATAL("The function index %i has not been found.", fct_idx);
                }
                expr += it->first;
            } else {
                expr += _named_functions[fct_idx].getName() + "(";
            }
            break;
        }
        case TraversePosition::BetweenChildren:
            expr += ", ";
            break;
        case TraversePosition::EndNode:
            expr += ")";
        }
    };

    traverse(_map_sink_source, fct_it->second, callback);
    DBUG("expression: %s", expr.c_str());
    return expr;
}

SpecificFunctionCaller
NamedFunctionCaller::getSpecificFunctionCaller(const std::string &function_name)
{
    auto const fct_it = _map_name_idx.find(function_name);
    if (fct_it == _map_name_idx.end()) {
        OGS_FATAL("A function with the name `%s' has not been found.",
                  function_name.c_str());
    }
    return SpecificFunctionCaller(fct_it->second, *this);
}

std::size_t NamedFunctionCaller::getNumberOfUnboundArguments() const
{
    return -_uninitialized - 1;
}

SpecificFunctionCaller::SpecificFunctionCaller(const std::size_t function_idx,
                                             const NamedFunctionCaller& caller)
    : _function_idx(function_idx), _caller(caller)
{
}

double SpecificFunctionCaller::call(
    const std::vector<double>& unbound_arguments) const
{
    return _caller.call(_function_idx, unbound_arguments);
}

std::size_t SpecificFunctionCaller::getNumberOfUnboundArguments() const
{
    return _caller.getNumberOfUnboundArguments();
}

} // namespace NumLib
