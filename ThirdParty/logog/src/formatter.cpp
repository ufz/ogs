 
/* 
 * \file formatter.cpp
 */

#include "logog.hpp"

namespace logog {

	Formatter::Formatter() :
		m_bShowTimeOfDay( false )
	{
		m_sMessageBuffer.reserve( LOGOG_FORMATTER_MAX_LENGTH );
		m_sIntBuffer.reserve_for_int();
	}

	void Formatter::RenderTimeOfDay()
	{
		if ( m_bShowTimeOfDay )
		{
#ifndef LOGOG_UNICODE
			TimeStamp stamp;
			m_sMessageBuffer.append( stamp.Get() );
			m_sMessageBuffer.append(": ");
#endif
		}
	}

	const LOGOG_CHAR * Formatter::ErrorDescription( const LOGOG_LEVEL_TYPE level )
	{
		if ( level <= LOGOG_LEVEL_NONE )
			return LOGOG_CONST_STRING("none");

		if ( level <= LOGOG_LEVEL_EMERGENCY )
			return LOGOG_CONST_STRING("emergency");

		if ( level <= LOGOG_LEVEL_ALERT )
			return LOGOG_CONST_STRING("alert");

		if ( level <= LOGOG_LEVEL_CRITICAL )
			return LOGOG_CONST_STRING("critical");

		if ( level <= LOGOG_LEVEL_ERROR )
			return LOGOG_CONST_STRING("error");

		if ( level <= LOGOG_LEVEL_WARN )
			return LOGOG_CONST_STRING("warning");

		if ( level <= LOGOG_LEVEL_INFO )
			return LOGOG_CONST_STRING("info");

		if ( level <= LOGOG_LEVEL_DEBUG )
			return LOGOG_CONST_STRING("debug");

		return LOGOG_CONST_STRING("unknown");
	}

	bool Formatter::GetShowTimeOfDay() const
	{
		return m_bShowTimeOfDay;
	}

	void Formatter::SetShowTimeOfDay( bool val )
	{
		m_bShowTimeOfDay = val;
	}

	TOPIC_FLAGS Formatter::GetTopicFlags( const Topic &topic )
	{
		return topic.GetTopicFlags();
	}

	LOGOG_STRING &FormatterGCC::Format( const Topic &topic, const Target &target )
	{
		TOPIC_FLAGS flags;
		flags = GetTopicFlags( topic );

		m_sMessageBuffer.clear();

		if ( flags & TOPIC_FILE_NAME_FLAG )
		{
			m_sMessageBuffer.append( topic.FileName() );
			m_sMessageBuffer.append( ':' );
		}

		if ( flags & TOPIC_LINE_NUMBER_FLAG )
		{
			m_sIntBuffer.assign( topic.LineNumber() );
			m_sMessageBuffer.append( m_sIntBuffer );

			m_sMessageBuffer.append( LOGOG_CONST_STRING(": "));
		}

		RenderTimeOfDay();

		if ( flags & TOPIC_LEVEL_FLAG )
		{
			m_sMessageBuffer.append( ErrorDescription( topic.Level()));
			m_sMessageBuffer.append( LOGOG_CONST_STRING(": "));
		}

		if ( flags & TOPIC_GROUP_FLAG )
		{
			m_sMessageBuffer.append( LOGOG_CONST_STRING("{") );
			m_sMessageBuffer.append( topic.Group() );
			m_sMessageBuffer.append( LOGOG_CONST_STRING("} ") );
		}

		if ( flags & TOPIC_CATEGORY_FLAG )
		{
			m_sMessageBuffer.append( LOGOG_CONST_STRING("["));
			m_sMessageBuffer.append( topic.Category() );
			m_sMessageBuffer.append( LOGOG_CONST_STRING("] "));
		}

		if ( flags & TOPIC_MESSAGE_FLAG )
		{
			m_sMessageBuffer.append( topic.Message() );
			m_sMessageBuffer.append( (LOGOG_CHAR)'\n' );
		}

		if ( target.GetNullTerminatesStrings() )
			m_sMessageBuffer.append( (LOGOG_CHAR)NULL );

		return m_sMessageBuffer;
	}


	LOGOG_STRING &FormatterMSVC::Format( const Topic &topic, const Target &target )
    {
        m_sMessageBuffer.clear();

        TOPIC_FLAGS flags;
        flags = GetTopicFlags( topic );

        if ( flags & TOPIC_FILE_NAME_FLAG )
        {
            m_sMessageBuffer.append( topic.FileName() );
            m_sMessageBuffer.append( '(' );
        }

        if ( flags & TOPIC_LINE_NUMBER_FLAG  )
        {
            m_sIntBuffer.assign( topic.LineNumber() );
            m_sMessageBuffer.append( m_sIntBuffer );

            m_sMessageBuffer.append( LOGOG_CONST_STRING(") : ") );
        }

		RenderTimeOfDay();

        if ( flags & TOPIC_LEVEL_FLAG )
        {
            m_sMessageBuffer.append( ErrorDescription( topic.Level() ) );
            m_sMessageBuffer.append( LOGOG_CONST_STRING(": "));
        }

        if ( flags & TOPIC_GROUP_FLAG )
        {
            m_sMessageBuffer.append( LOGOG_CONST_STRING("{"));
            m_sMessageBuffer.append( topic.Group() );
            m_sMessageBuffer.append( LOGOG_CONST_STRING("} "));
        }

        if ( flags & TOPIC_CATEGORY_FLAG )
        {
            m_sMessageBuffer.append( LOGOG_CONST_STRING("["));
            m_sMessageBuffer.append( topic.Category() );
            m_sMessageBuffer.append( LOGOG_CONST_STRING("] "));
        }

        if ( flags & TOPIC_MESSAGE_FLAG )
        {
            m_sMessageBuffer.append( topic.Message() );
#ifdef LOGOG_FLAVOR_WINDOWS
			m_sMessageBuffer.append( LOGOG_CONST_STRING("\r\n") );
#else // LOGOG_FLAVOR_WINDOWS
            m_sMessageBuffer.append( LOGOG_CONST_STRING("\n") );
#endif // LOGOG_FLAVOR_WINDOWS
        }

		if ( target.GetNullTerminatesStrings() )
			m_sMessageBuffer.append( LOGOG_CHAR( NULL ) );

        return m_sMessageBuffer;
    }

	Formatter &GetDefaultFormatter()
	{
		Statics *pStatic = &Static();

		if ( pStatic->s_pDefaultFormatter == NULL )
		{
#ifdef LOGOG_FLAVOR_WINDOWS
			pStatic->s_pDefaultFormatter = new FormatterMSVC();
#else
			pStatic->s_pDefaultFormatter = new FormatterGCC();
#endif
		}

		return *( pStatic->s_pDefaultFormatter );
	}

	void DestroyDefaultFormatter()
	{
		Statics *pStatic = &Static();
		Formatter *pDefaultFormatter = pStatic->s_pDefaultFormatter;

		if ( pDefaultFormatter != NULL )
			delete pDefaultFormatter;

		pStatic->s_pDefaultFormatter = NULL;
}

const char * TimeStamp::Get()
{
	time_t tRawTime;
	struct tm * tmInfo;

	time ( &tRawTime );

#ifdef LOGOG_FLAVOR_WINDOWS
#pragma warning( push )
#pragma warning( disable : 4996 )
#endif // LOGOG_FLAVOR_WINDOWS
	/* Microsoft is afraid of this function; I'm not sure this warning is sensible */
	tmInfo = localtime ( &tRawTime );
#ifdef LOGOG_FLAVOR_WINDOWS
#pragma warning( pop )
#endif // LOGOG_FLAVOR_WINDOWS

	strftime (cTimeString, LOGOG_TIME_STRING_MAX, "%c", tmInfo);

	return cTimeString;
}

}

