/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#pragma once

#include <string>
#include <vector>
#include <list>
#include <iosfwd>
#include <algorithm>

#include <tclap/CmdLineInterface.h>
#include <tclap/CmdLineOutput.h>
#include <tclap/StdOutput.h>
#include <tclap/XorHandler.h>
#include <tclap/Arg.h>

namespace BaseLib
{

/**
 * TCLAP standard output modified as follows
 * - Print arguments in the order of added to Command object
 */
class TCLAPCustomOutput : public TCLAP::StdOutput
{
public:
    /**
     * Prints the usage to stdout.  Can be overridden to
     * produce alternative behavior.
     * \param c - The CmdLine object the output is generated for.
     */
    virtual void usage(TCLAP::CmdLineInterface& c);

    /**
     * Prints (to stderr) an error message, short usage
     * Can be overridden to produce alternative behavior.
     * \param c - The CmdLine object the output is generated for.
     * \param e - The ArgException that caused the failure.
     */
    virtual void failure(TCLAP::CmdLineInterface& c,
            TCLAP::ArgException& e );

protected:

    /**
     * Writes a brief usage message with short args.
     * \param c - The CmdLine object the output is generated for.
     * \param os - The stream to write the message to.
     */
    void _shortUsage( TCLAP::CmdLineInterface& c, std::ostream& os ) const;

    /**
     * Writes a longer usage message with long and short args,
     * provides descriptions and prints message.
     * \param c - The CmdLine object the output is generated for.
     * \param os - The stream to write the message to.
     */
    void _longUsage( TCLAP::CmdLineInterface& c, std::ostream& os ) const;

};

inline void TCLAPCustomOutput::usage(TCLAP::CmdLineInterface& _cmd )
{
    std::cout << std::endl << "USAGE: " << std::endl << std::endl;

    _shortUsage( _cmd, std::cout );

    std::cout << std::endl << std::endl << "Where: " << std::endl << std::endl;

    _longUsage( _cmd, std::cout );

    std::cout << std::endl;

}

inline void TCLAPCustomOutput::failure( TCLAP::CmdLineInterface& _cmd,
        TCLAP::ArgException& e )
{
    std::string progName = _cmd.getProgramName();

    std::cerr << "PARSE ERROR: " << e.argId() << std::endl
              << "             " << e.error() << std::endl << std::endl;

    if ( _cmd.hasHelpAndVersion() )
        {
            std::cerr << "Brief USAGE: " << std::endl;

            _shortUsage( _cmd, std::cerr );

            std::cerr << std::endl << "For complete USAGE and HELP type: "
                      << std::endl << "   " << progName << " --help"
                      << std::endl << std::endl;
        }
    else
        usage(_cmd);

    throw TCLAP::ExitException(1);
}

inline void
TCLAPCustomOutput::_shortUsage( TCLAP::CmdLineInterface& _cmd,
                        std::ostream& os ) const
{
    std::list<TCLAP::Arg*> argList = _cmd.getArgList();
    std::string progName = _cmd.getProgramName();
    TCLAP::XorHandler xorHandler = _cmd.getXorHandler();
    std::vector< std::vector<TCLAP::Arg*> > xorList = xorHandler.getXorList();

    std::string s = progName + " ";

    // first the xor
    for ( int i = 0; static_cast<unsigned int>(i) < xorList.size(); i++ )
        {
            s += " {";
            for ( TCLAP::ArgVectorIterator it = xorList[i].begin();
                  it != xorList[i].end(); it++ )
                s += (*it)->shortID() + "|";

            s[s.length()-1] = '}';
        }

    // then the rest
    for (auto it = argList.rbegin(); it != argList.rend(); ++it) // here modified
        if ( !xorHandler.contains( (*it) ) )
            s += " " + (*it)->shortID();

    // if the program name is too long, then adjust the second line offset
    int secondLineOffset = static_cast<int>(progName.length()) + 2;
    if ( secondLineOffset > 75/2 )
        secondLineOffset = static_cast<int>(75/2);

    spacePrint( os, s, 75, 3, secondLineOffset );
}

inline void
TCLAPCustomOutput::_longUsage( TCLAP::CmdLineInterface& _cmd,
                       std::ostream& os ) const
{
    std::list<TCLAP::Arg*> argList = _cmd.getArgList();
    std::string message = _cmd.getMessage();
    TCLAP::XorHandler xorHandler = _cmd.getXorHandler();
    std::vector< std::vector<TCLAP::Arg*> > xorList = xorHandler.getXorList();

    // first the xor
    for ( int i = 0; static_cast<unsigned int>(i) < xorList.size(); i++ )
        {
            for ( TCLAP::ArgVectorIterator it = xorList[i].begin();
                  it != xorList[i].end();
                  it++ )
                {
                    this->spacePrint( os, (*it)->longID(), 75, 3, 3 );
                    spacePrint( os, (*it)->getDescription(), 75, 5, 0 );

                    if ( it+1 != xorList[i].end() )
                        spacePrint(os, "-- OR --", 75, 9, 0);
                }
            os << std::endl << std::endl;
        }

    // then the rest
    for (auto it = argList.rbegin(); it != argList.rend(); it++) // here modified
        if ( !xorHandler.contains( (*it) ) )
            {
                spacePrint( os, (*it)->longID(), 75, 3, 3 );
                spacePrint( os, (*it)->getDescription(), 75, 5, 0 );
                os << std::endl;
            }

    os << std::endl;

    spacePrint( os, message, 75, 3, 0 );
}

} //namespace BaseLib
