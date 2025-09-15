/**
 * \file
 * \copyright
 * Copyright (c) 2025-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include <tclap/CmdLine.h>

#include <regex>

namespace BaseLib
{

template <typename T>
const char* type_name();

template <>
const char* type_name<std::string>()
{
    return "string";
}
template <>
const char* type_name<int>()
{
    return "int";
}
template <>
const char* type_name<double>()
{
    return "double";
}
template <>
const char* type_name<float>()
{
    return "float";
}
template <>
const char* type_name<bool>()
{
    return "bool";
}
template <>
const char* type_name<size_t>()
{
    return "size_t";
}
template <>
const char* type_name<long>()
{
    return "long";
}
template <>
const char* type_name<unsigned>()
{
    return "unsigned";
}

class TCLAPOutput : public TCLAP::StdOutput
{
public:
    virtual void usage(TCLAP::CmdLineInterface& _cmd)
    {
        // Follows implementation from tclap
        std::cout << std::endl << "USAGE: " << std::endl << std::endl;
        _shortUsage(_cmd, std::cout);
        std::cout << std::endl
                  << std::endl
                  << "Where: " << std::endl
                  << std::endl;

        std::string message = _cmd.getMessage();
        std::list<TCLAP::ArgGroup*> argSets = _cmd.getArgGroups();

        std::list<TCLAP::Arg*> unlabelled;
        for (std::list<TCLAP::ArgGroup*>::iterator sit = argSets.begin();
             sit != argSets.end();
             ++sit)
        {
            TCLAP::ArgGroup& argGroup = **sit;

            int visible = CountVisibleArgs(argGroup);
            bool exclusive = visible > 1 && argGroup.isExclusive();
            bool forceRequired = visible == 1 && argGroup.isRequired();
            if (exclusive)
            {
                spacePrint(std::cout,
                           argGroup.isRequired() ? "One of:" : "Either of:", 75,
                           3, 0);
            }

            for (TCLAP::ArgGroup::iterator it = argGroup.begin();
                 it != argGroup.end();
                 ++it)
            {
                TCLAP::Arg& arg = **it;

                if (!arg.visibleInHelp())
                {
                    continue;
                }

                if (!arg.hasLabel())
                {
                    unlabelled.push_back(&arg);
                    continue;
                }

                // Added by LB:
                std::string shortID = arg.shortID();
                std::smatch match;
                bool required = arg.isRequired() || forceRequired;
                std::stringstream ss;
                if (std::regex_search(shortID, match, std::regex("<([^>]+)>")))
                {
                    std::string identifier = match[1];

                    ss << identifier;
                    tryPrintValue<std::string>(*it, ss, required);
                    tryPrintValue<int>(*it, ss, required);
                    tryPrintValue<double>(*it, ss, required);
                    tryPrintValue<float>(*it, ss, required);
                    tryPrintValue<size_t>(*it, ss, required);
                    tryPrintValue<unsigned>(*it, ss, required);
                    tryPrintValue<unsigned int>(*it, ss, required);
                    ss << ": ";
                }
                ss << arg.getDescription(required);

                if (exclusive)
                {
                    spacePrint(std::cout, arg.longID(), 75, 6, 3);
                    spacePrint(std::cout, ss.str(), 75, 8, 0);
                }
                else
                {
                    spacePrint(std::cout, arg.longID(), 75, 3, 3);
                    spacePrint(std::cout, ss.str(), 75, 5, 0);
                }
                std::cout << '\n';
            }
        }

        for (TCLAP::ArgListIterator it = unlabelled.begin();
             it != unlabelled.end();
             ++it)
        {
            const TCLAP::Arg& arg = **it;
            spacePrint(std::cout, arg.longID(), 75, 3, 3);
            spacePrint(std::cout, arg.getDescription(), 75, 5, 0);
            std::cout << '\n';
        }

        if (!message.empty())
        {
            spacePrint(std::cout, message, 75, 3, 0);
        }

        std::cout.flush();
    }

private:
    template <typename T>
    void tryPrintValue(TCLAP::Arg* arg, std::ostream& os, bool required)
    {
        if (auto val = dynamic_cast<TCLAP::ValueArg<T>*>(arg))
        {
            os << " (" << type_name<T>();
            if (!required)
            {
                if constexpr (std::is_same_v<T, std::string>)
                {
                    // print non-empty default value only
                    if (val->getValue() != "")
                    {
                        os << ", default: " << val->getValue();
                    }
                }
                else
                {
                    os << ", default: " << val->getValue();
                }
            }
            os << ")";
        }
    }
};

}  // namespace BaseLib
