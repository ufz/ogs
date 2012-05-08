 
/* 
 * \file topic.cpp
 */

#include "logog.hpp"

namespace logog {

	void SetDefaultLevel( LOGOG_LEVEL_TYPE level )
	{
		Filter *pDefaultFilter = &GetDefaultFilter();

		pDefaultFilter->Level( level );
	}

	Topic::Topic( const LOGOG_LEVEL_TYPE level ,
		const LOGOG_CHAR *sFileName ,
		const int nLineNumber ,
		const LOGOG_CHAR *sGroup ,
		const LOGOG_CHAR *sCategory ,
		const LOGOG_CHAR *sMessage ,
		const double dTimestamp )
	{
		m_TopicFlags = 0;

		if ( sFileName != NULL )
		{
			m_vStringProps[ TOPIC_FILE_NAME ] = sFileName;
			m_TopicFlags |= TOPIC_FILE_NAME_FLAG;
		}

		if ( sGroup != NULL )
		{
			m_vStringProps[ TOPIC_GROUP ] = sGroup;
			m_TopicFlags |= TOPIC_GROUP_FLAG;
		}

		if ( sCategory != NULL )
		{
			m_vStringProps[ TOPIC_CATEGORY ] = sCategory;
			m_TopicFlags |= TOPIC_CATEGORY_FLAG;
		}

		if ( sMessage != NULL )
		{
			m_vStringProps[ TOPIC_MESSAGE ] = sMessage;
			m_TopicFlags |= TOPIC_MESSAGE_FLAG;
		}

		m_vIntProps[ TOPIC_LEVEL ] = level;

		if ( level != LOGOG_LEVEL_ALL )
		{
			m_TopicFlags |= TOPIC_LEVEL_FLAG;
		}

		m_vIntProps[ TOPIC_LINE_NUMBER ] = nLineNumber;

		if ( nLineNumber != 0 )
		{
			m_TopicFlags |= TOPIC_LINE_NUMBER_FLAG;
		}

		m_tTime = dTimestamp;

		if ( dTimestamp != 0.0f ) //-V550
			m_TopicFlags |= TOPIC_TIMESTAMP_FLAG;
	}

	bool Topic::IsTopic() const
	{
		return true;
	}

	int Topic::Send( const Topic &node )
	{
		LockableNodesType::iterator it;

		{
			ScopedLock sl( m_Subscribers );
			it = m_Subscribers.begin();
		}

		/* Iterate over the subscribers, but only addressing the subscribers group while locking it */
		Topic *pCurrentTopic;
		Node *pCurrentNode;
		m_Subscribers.MutexLock();
		int nError = 0;

		while ( it != m_Subscribers.end() )
		{
			pCurrentNode = *it;

			if ( pCurrentNode->IsTopic() == false )
				continue;

			pCurrentTopic = ( Topic * )pCurrentNode;

			if ( pCurrentTopic )
				nError += pCurrentTopic->Receive( node );

			it++;
		}

		m_Subscribers.MutexUnlock();

		return nError;
	}

	int Topic::Transmit()
	{
		return Send( *this );
	}

	int Topic::Receive( const Topic &node )
	{
		/* Default implementation -- send it on to all children */
		return Send( node );
	}

	bool Topic::CanSubscribeTo( const Node &otherNode )
	{
		if ( CanSubscribe() == false )
			return false;

		if ( otherNode.IsTopic() == false )
			return false;

		Topic *pTopic = ( Topic * )&otherNode;

		/* This function will change from topic class to topic class. */
		return CanSubscribeCheckTopic( *pTopic );
	}

	bool Topic::CanSubscribeCheckTopic( const Topic &other )
	{
		/* This is the generic comparison case.  We'll want to optimize this function for other types
		* of topics.
		*/

		/* Check topics in likely order of disinterest */
		if ( m_TopicFlags & TOPIC_LEVEL_FLAG )
		{
			/* Topic levels are less interesting the larger the numbers are. */
			if ( other.m_vIntProps[ TOPIC_LEVEL ] > m_vIntProps[ TOPIC_LEVEL ] )
				return false;
		}

		if ( m_TopicFlags & TOPIC_GROUP_FLAG )
		{
			/* If our topic is not a substring of the publisher's topic, ignore this */
			if (( other.m_vStringProps[ TOPIC_GROUP ] ).find( m_vStringProps[ TOPIC_GROUP ] ) == LOGOG_STRING::npos )
				return false;
		}

		if ( m_TopicFlags & TOPIC_CATEGORY_FLAG )
		{
			/* If our topic is not a substring of the publisher's topic, ignore this */
			if (( other.m_vStringProps[ TOPIC_CATEGORY ] ).find( m_vStringProps[ TOPIC_CATEGORY ] ) == LOGOG_STRING::npos )
				return false;
		}

		if ( m_TopicFlags & TOPIC_FILE_NAME_FLAG )
		{
			/* If our topic is not a substring of the publisher's file name, ignore this. */
			if (( other.m_vStringProps[ TOPIC_FILE_NAME ] ).find( m_vStringProps[ TOPIC_FILE_NAME ] ) == LOGOG_STRING::npos )
				return false;
		}

		if ( m_TopicFlags & TOPIC_LINE_NUMBER_FLAG )
		{
			/* If our line number doesn't equal theirs, ignore this */
			if ( other.m_vIntProps[ TOPIC_LINE_NUMBER ] != m_vIntProps[ TOPIC_LINE_NUMBER ] )
				return false;
		}

		if ( m_TopicFlags & TOPIC_MESSAGE_FLAG )
		{
			/* If our topic is not a substring of the publisher's file name, ignore this. */
			if (( other.m_vStringProps[ TOPIC_MESSAGE ] ).find( m_vStringProps[ TOPIC_MESSAGE ] ) == LOGOG_STRING::npos )
				return false;
		}

		if ( m_TopicFlags & TOPIC_TIMESTAMP_FLAG )
		{
			/* Timestamps are only interesting if they're greater than or equal to ours. */
			if ( other.m_tTime < m_tTime )
				return false;
		}

		/* all tests passed */
		return true;
	}

	bool Topic::PublishTo( Node &subscriber )
	{
#ifdef LOGOG_INTERNAL_DEBUGGING
		if ( &subscriber == this )
			LOGOG_INTERNAL_FAILURE;
#endif
		bool bWasInserted;

		/** Additional checking may be required first -- can the subscriber handle this publishing? */
		if ( subscriber.IsTopic() )
		{
			Topic *pSubscriber = (Topic *)&subscriber;

			if ( pSubscriber->CanSubscribeTo( *this ) == false )
				return false;
		}

		{
			ScopedLock sl( m_Subscribers );
			bWasInserted = ( m_Subscribers.insert( &subscriber ) ).second;
		}

		if ( bWasInserted )
			subscriber.SubscribeTo( *this );

		return bWasInserted;
	}

	void Topic::Format( const LOGOG_CHAR *cFormatMessage, ... )
	{
		va_list args;

		va_start( args, cFormatMessage );
		m_vStringProps[ TOPIC_MESSAGE ].format_va( cFormatMessage, args );
		va_end( args );

		m_TopicFlags |= TOPIC_MESSAGE_FLAG;
	}

	const LOGOG_STRING & Topic::FileName() const
	{
		return m_vStringProps[ TOPIC_FILE_NAME ];
	}

	void Topic::FileName( const LOGOG_STRING &s )
	{
		m_vStringProps[ TOPIC_FILE_NAME ] = s;
		m_TopicFlags |= TOPIC_FILE_NAME_FLAG;
	}

	const LOGOG_STRING & Topic::Message() const
	{
		return m_vStringProps[ TOPIC_MESSAGE ];
	}

	void Topic::Message( const LOGOG_STRING &s )
	{
		m_vStringProps[ TOPIC_MESSAGE ] = s;
		m_TopicFlags |= TOPIC_MESSAGE_FLAG;
	}

	const LOGOG_STRING & Topic::Category() const
	{
		return m_vStringProps[ TOPIC_CATEGORY ];
	}

	void Topic::Category( const LOGOG_STRING &s )
	{
		m_vStringProps[ TOPIC_CATEGORY ] = s;
		m_TopicFlags |= TOPIC_CATEGORY_FLAG;
	}

	const LOGOG_STRING & Topic::Group() const
	{
		return m_vStringProps[ TOPIC_GROUP ];
	}

	void Topic::Group( const LOGOG_STRING &s )
	{
		m_vStringProps[ TOPIC_GROUP ] = s;
		m_TopicFlags |= TOPIC_GROUP_FLAG;
	}

	int Topic::LineNumber() const
	{
		return m_vIntProps[ TOPIC_LINE_NUMBER ];
	}

	void Topic::LineNumber( const int num )
	{
		m_vIntProps[ TOPIC_LINE_NUMBER ] = num;
		m_TopicFlags |= TOPIC_LINE_NUMBER_FLAG;
	}

	LOGOG_LEVEL_TYPE Topic::Level() const
	{
		return ( LOGOG_LEVEL_TYPE )m_vIntProps[ TOPIC_LEVEL ];
	}

	void Topic::Level( LOGOG_LEVEL_TYPE level )
	{
		m_vIntProps[ TOPIC_LEVEL ] = level;
		m_TopicFlags |= TOPIC_LEVEL_FLAG;
	}

	logog::LOGOG_TIME Topic::Timestamp() const
	{
		return m_tTime;
	}

	void Topic::Timestamp( const LOGOG_TIME t )
	{
		m_tTime = t;
		m_TopicFlags |= TOPIC_TIMESTAMP_FLAG;
	}

	TOPIC_FLAGS Topic::GetTopicFlags() const
	{
		return m_TopicFlags;
	}


/********************************************************/

	Filter::Filter( const LOGOG_LEVEL_TYPE level ,
		const LOGOG_CHAR *sFileName ,
		const int nLineNumber ,
		const LOGOG_CHAR *sGroup ,
		const LOGOG_CHAR *sCategory ,
		const LOGOG_CHAR *sMessage ,
		const double dTimestamp ) :
	Topic( level, sFileName, nLineNumber, sGroup, sCategory, sMessage, dTimestamp )
	{
		Statics *pStatic = &Static();

#ifdef LOGOG_INTERNAL_DEBUGGING
		if ( pStatic == NULL )
			LOGOG_INTERNAL_FAILURE;
#endif

		if ( pStatic->s_pDefaultFilter == NULL )
			pStatic->s_pDefaultFilter = this;

		PublishToMultiple( AllTargets() );

		LockableNodesType *pFilterNodes = &AllFilters();

		{
			ScopedLock sl( *pFilterNodes );
			pFilterNodes->insert( this );
		}
	}

	Filter &GetDefaultFilter()
	{
		Statics *pStatic = &Static();

#ifdef LOGOG_INTERNAL_DEBUGGING
		if ( pStatic == NULL )
			LOGOG_INTERNAL_FAILURE;
#endif

		if ( pStatic->s_pDefaultFilter == NULL )
		{
			pStatic->s_pDefaultFilter = new Filter( LOGOG_LEVEL );
		}

		return *((Filter *)(pStatic->s_pDefaultFilter));
	}

	TopicGroup::TopicGroup( const LOGOG_CHAR *sGroup ) :
		Topic( LOGOG_LEVEL_ALL, NULL, 0, sGroup )
	{
	}

	bool TopicGroup::CanSubscribeCheckTopic( const Topic &other )
	{
		if ( m_TopicFlags & TOPIC_LEVEL_FLAG )
		{
			/* Topic levels are less interesting the larger the numbers are. */
			if ( other.m_vIntProps[ TOPIC_LEVEL ] > m_vIntProps[ TOPIC_LEVEL ] )
				return false;
		}

		return true;
	}

	TopicLevel::TopicLevel( const LOGOG_LEVEL_TYPE level ) :
	Topic( level )
	{
	}


	bool TopicLevel::CanSubscribeCheckTopic( const Topic &other )
	{
		/* Check topics in likely order of disinterest */
		if ( m_TopicFlags & TOPIC_LEVEL_FLAG )
		{
			/* Topic levels are less interesting the larger the numbers are. */
			if ( other.m_vIntProps[ TOPIC_LEVEL ] > m_vIntProps[ TOPIC_LEVEL ] )
				return false;
		}

		/* all tests passed */
		return true;
	}

	TopicSource::TopicSource( const LOGOG_LEVEL_TYPE level ,
		const LOGOG_CHAR *sFileName,
		const int nLineNumber,
		const LOGOG_CHAR *sGroup,
		const LOGOG_CHAR *sCategory,
		const LOGOG_CHAR *sMessage,
		const double dTimestamp ) :
	Topic( level, sFileName, nLineNumber, sGroup, sCategory, sMessage, dTimestamp )
	{
	}

	bool TopicSource::SubscribeTo( Node & )
	{
		return false;
	}

	bool TopicSource::UnsubscribeTo( Node & )
	{
		return false;
	}

	bool TopicSource::CanSubscribe() const
	{
		return false;
	}

	bool TopicSink::IsTopic() const
	{
		return true;
	}

	void TopicSink::Initialize()
	{

	}

	bool TopicSink::PublishTo( Node & )
	{
		return false;
	}

	bool TopicSink::UnpublishTo( Node & )
	{
		return false;
	}

	bool TopicSink::CanPublish() const
	{
		return false;
	}
}

