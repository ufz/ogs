/**
 * \file checkpoint.hpp Representations of a program counter reaching a specific point in code.
 */

#ifndef __LOGOG_CHECKPOINT_HPP__
#define __LOGOG_CHECKPOINT_HPP__

namespace logog
{
/** A checkpoint is a topic that fires when a specific section of code is executed.  The first time a bit of
 ** code is executed, a Checkpoint is instanced, and when the code is executed again, the Checkpoint is
 ** reused.
 **/
class Checkpoint : public TopicSource
{
public:
    Checkpoint( const LOGOG_LEVEL_TYPE level = LOGOG_LEVEL_ALL,
                const LOGOG_CHAR *sFileName = NULL,
                const int nLineNumber = 0,
                const LOGOG_CHAR *sGroup = NULL,
                const LOGOG_CHAR *sCategory = NULL,
                const LOGOG_CHAR *sMessage = NULL,
                const double dTimestamp = 0.0f );

    /** Sends the node in question.  Optionally updates the timestamp in this checkpoint before sending the node. */
    virtual int Send( const Topic &node );

};
}

#endif // __LOGOG_CHECKPOINT_HPP_
