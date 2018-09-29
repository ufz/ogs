/**
 * \file message.hpp Messages; items transmitted to a log.
 */

#ifndef __LOGOG_MESSAGE_HPP__
#define __LOGOG_MESSAGE_HPP__

namespace logog
{

/** A message is a piece of text that's actually transmitted to outputs.  Messages can be asked to 
 ** Transmit() themselves once they are created.
 **/
class Message : public Checkpoint
{
public:
    Message( const LOGOG_LEVEL_TYPE level = LOGOG_LEVEL_ALL,
             const LOGOG_CHAR *sFileName = NULL,
             const int nLineNumber = 0,
             const LOGOG_CHAR *sGroup = NULL,
             const LOGOG_CHAR *sCategory = NULL,
             const LOGOG_CHAR *sMessage = NULL,
             const double dTimestamp = 0.0f,
			 bool *bIsCreated = NULL );

	virtual ~Message();

    /** Causes this checkpoint to republish itself to all existing filters after
      * unpublishing itself.  This can be necessary if the message within this
      * message has changed in such a way that the downstream Filter objects
      * might react differently to it.
      */
    virtual bool Republish();

	Mutex m_Transmitting;
	bool *m_pbIsCreated;
};

extern Mutex &GetMessageCreationMutex();
extern void DestroyMessageCreationMutex();

}


#endif // __LOGOG_MESSAGE_HPP_
