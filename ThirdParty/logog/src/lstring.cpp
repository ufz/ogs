 /* 
 * \file lstring.cpp
 */

#include "logog.hpp"

namespace logog {

	String::String()
	{
		Initialize();
	}

	String::~String()
	{
		Free();
	}

	void String::Free()
	{
		if ( m_pBuffer && ( m_bIsConst == false ))
		{
			Deallocate( (void *)m_pBuffer );
			m_pBuffer = m_pEndOfBuffer = m_pOffset = NULL;
		}

		if ( m_pKMP )
		{
			Deallocate( (void *)m_pKMP );
			m_pKMP = NULL;
		}
	}

	size_t String::Length( const LOGOG_CHAR *chars )
	{
		unsigned int len = 0;

		while ( *chars++ )
			len++;

		return len;
	}

	String & String::operator=( const String & other )
	{
		Free();
		Initialize();
		assign( other );
		return *this;
	}

	String & String::operator=( const LOGOG_CHAR *pstr )
	{
		Free();
		Initialize();
		assign( pstr );
		return *this;
	}

	String::String( const String &other )
	{
		Initialize();
		assign( other );
	}

	String::String( const LOGOG_CHAR *pstr )
	{
		Initialize();
		assign( pstr );
	}

	size_t String::size() const
	{
		return ( m_pOffset - m_pBuffer );
	}

	void String::clear()
	{
		m_pOffset = m_pBuffer;
	}

	size_t String::reserve( size_t nSize )
	{
		if ( nSize == (unsigned int)( m_pOffset - m_pBuffer ))
			return nSize;

		if ( nSize == 0 )
		{
			if ( m_pBuffer && ( *m_pBuffer != (LOGOG_CHAR)NULL ))
				Deallocate( (void *)m_pBuffer );

			Initialize();
			return 0;
		}

		LOGOG_CHAR *pNewBuffer = (LOGOG_CHAR *)Allocate( sizeof( LOGOG_CHAR ) * nSize );
		LOGOG_CHAR *pNewOffset = pNewBuffer;
		LOGOG_CHAR *pNewEnd = pNewBuffer + nSize;

		LOGOG_CHAR *pOldOffset = m_pOffset;

		if ( pOldOffset != NULL )
		{
			while (( pNewOffset < pNewEnd ) && ( *pOldOffset != (LOGOG_CHAR)NULL ))
				*pNewOffset++ = *pOldOffset++;
		}

		if (( m_pBuffer != NULL ) && ( m_bIsConst == false ))
			Deallocate( m_pBuffer );

		m_pBuffer = pNewBuffer;
		m_pOffset = pNewBuffer;
		m_pEndOfBuffer = pNewEnd;

		return ( m_pOffset - m_pBuffer );
	}

	size_t String::reserve_for_int()
	{
		reserve( 32 );
		return 32;
	}

	String::operator const LOGOG_CHAR *() const
	{
		return m_pBuffer;
	}

	const LOGOG_CHAR* String::c_str() const
	{
		return m_pBuffer;
	}

	size_t String::assign( const String &other )
	{
#ifdef LOGOG_INTERNAL_DEBUGGING
		if ( m_bIsConst )
			cout << "Can't reassign const string!" << endl;
#endif
		LOGOG_CHAR *pOther = other.m_pBuffer;

		if ( pOther == NULL )
			return 0;

		size_t othersize = other.size();

		reserve( othersize + 1 );
		m_pOffset = m_pBuffer;

		for ( unsigned int t = 0; t <= othersize ; t++ )
			*m_pOffset++ = *pOther++;

		return this->size();
	}

	size_t String::assign( const int value )
	{
#ifdef LOGOG_INTERNAL_DEBUGGING
		if ( m_bIsConst )
			cout << "Can't reassign const string!" << endl;
#endif

		int number = value;
		m_pOffset = m_pBuffer;

		int bSign = value;

		if (( bSign = number) < 0)
			number = -number;

		do
		{
			*m_pOffset++ = _LG("0123456789")[ number % 10 ];
		}
		while( number /= 10 );

		if (bSign < 0)
			*m_pOffset++ = '-';

		*m_pOffset = (LOGOG_CHAR)'\0';

		reverse( m_pBuffer, m_pOffset - 1 );

		return ( m_pOffset - m_pBuffer );
	}

	size_t String::assign( const LOGOG_CHAR *other, const LOGOG_CHAR *pEnd /*= NULL */ )
	{
		size_t len;

		if ( pEnd == NULL )
			len = Length( other );
		else
			len = ( pEnd - other );
		/** This constant decides whether assigning a LOGOG_CHAR * to a String will cause the String to use the previous buffer
		* in place, or create a new buffer and copy the results.
		*/
#ifdef LOGOG_COPY_CONST_CHAR_ARRAY_ON_ASSIGNMENT
		reserve( len + 1 );

		for (unsigned int t = 0; t <= len; t++ )
			*m_pOffset++ = *other++;
#else  // LOGOG_COPY_CONST_CHAR_ARRAY_ON_ASSIGNMENT

#ifdef LOGOG_INTERNAL_DEBUGGING
		if ( m_bIsConst )
			cout << "Can't reassign const string!" << endl;
#endif
		/* In this case we don't copy the buffer, just reuse it */
		m_pBuffer = const_cast< LOGOG_CHAR *>( other );
		m_pOffset = m_pBuffer + len + 1;
		m_pEndOfBuffer = m_pOffset;
		m_bIsConst = true;

#endif // LOGOG_COPY_CONST_CHAR_ARRAY_ON_ASSIGNMENT

		return (int) len;
	}
	size_t String::append( const String &other )
	{
#ifdef LOGOG_INTERNAL_DEBUGGING
		if ( m_bIsConst )
			cout << "Can't reassign const string!" << endl;
#endif

		return append( other.m_pBuffer );
	}

	size_t String::append( const LOGOG_CHAR *other )
	{
#ifdef LOGOG_INTERNAL_DEBUGGING
		if ( m_bIsConst )
			cout << "Can't reassign const string!" << endl;
#endif
		if ( other == NULL )
			return 0;

		while (( m_pOffset < m_pEndOfBuffer ) && ( *other != (LOGOG_CHAR)NULL ))
			*m_pOffset++ = *other++;

		return ( m_pOffset - m_pBuffer );
	}

	size_t String::append( const LOGOG_CHAR c )
	{
		if ( m_pOffset < m_pEndOfBuffer )
			*m_pOffset++ = c;

		return ( m_pOffset - m_pBuffer );
	}
	void String::reverse( LOGOG_CHAR* pStart, LOGOG_CHAR* pEnd )
	{
		LOGOG_CHAR temp;

		while( pEnd > pStart)
		{
			temp=*pEnd, *pEnd-- =*pStart, *pStart++=temp;
		}
	}

	bool String::is_valid() const
	{
		return ( m_pBuffer != NULL );
	}

	size_t String::find( String &other ) const
	{
		if ( is_valid() && other.is_valid())
		{
			// KMP solution
			// String *pThis = const_cast< String *>(this);
			// return pThis->KMP( other.m_pBuffer, other.size());
			LOGOG_CHAR *pFound;

#ifdef LOGOG_UNICODE
			pFound = wcsstr( m_pBuffer, other.m_pBuffer );
#else // LOGOG_UNICODE
			pFound = strstr( m_pBuffer, other.m_pBuffer );
#endif

			if ( pFound != NULL )
			{
				return ( pFound - m_pBuffer );
			}

			return npos;
		}

		return npos;
	}

	void String::format( const LOGOG_CHAR *cFormatString, ... )
	{
		va_list args;

		va_start(args, cFormatString);
		format_va( cFormatString, args );
		va_end( args );
	}

	void String::format_va( const LOGOG_CHAR *cFormatString, va_list args )
	{
		int nActualSize = -1, nAttemptedSize;
		LOGOG_CHAR *pszFormatted = NULL;

		Free();

		/* Estimate length of output; don't pull in strlen() if we can help it */
		int nEstLength = 0;
		const LOGOG_CHAR *pCurChar = cFormatString;
		while ( *pCurChar++ )
			nEstLength++;

		if ( nEstLength == 0 )
		{
			clear();
			return;
		}

		/** nAttemptedSize is now a guess at an appropriate size, which is about 
		 ** two times the number of LOGOG_CHARs in the incoming format string.
		 **/
		nAttemptedSize = nEstLength * 2 * sizeof( LOGOG_CHAR );

		/* Some *printf implementations, such as msvc's, return -1 on failure.  
		 * Others, such as gcc, return the number
		 * of characters actually formatted on failure.  Deal with either case here.
		 */
		for ( ; ; )
		{
			/** We'll allocate that number of bytes.  NOTE that this has less of a chance
			 ** of working on a Unicode build.
			 **/
			pszFormatted = (LOGOG_CHAR *)Allocate( nAttemptedSize );
			if ( !pszFormatted )
			{
				LOGOG_INTERNAL_FAILURE;
			}

			*pszFormatted = (LOGOG_CHAR)'\0';

			va_list argsCopy;

			/** The va_list structure is not standardized across all platforms; in particular
			 ** Microsoft seems to have problem with the concept.
			 **/
#if defined( va_copy )
			va_copy( argsCopy, args );
#elif defined( __va_copy )
			__va_copy( argsCopy, args );
#else
			memcpy( &argsCopy, &args, sizeof(va_list) );
#endif

#ifdef LOGOG_UNICODE
			/** At this point, nSizeInWords will contain the number of words permitted in the
			 ** output buffer.  It takes into account space for appending a null character in the output
			 ** buffer as well.
			 **/
			int nSizeInWords = (nAttemptedSize / sizeof( LOGOG_CHAR ));
#endif
			/** The nActualSize value receives different things on different platforms.
			 ** On some platforms it receives -1 on failure; on other platforms
			 ** it receives the number of LOGOG_CHARs actually formatted (excluding
			 ** the trailing NULL).
			 **/

#ifdef LOGOG_FLAVOR_WINDOWS
#ifdef LOGOG_UNICODE
			nActualSize = _vsnwprintf_s( pszFormatted, nSizeInWords, _TRUNCATE, cFormatString, argsCopy );
#else // LOGOG_UNICODE
			nActualSize = vsnprintf_s( pszFormatted, nAttemptedSize, _TRUNCATE, cFormatString, argsCopy );
#endif // LOGOG_UNICODE
#else // LOGOG_FLAVOR_WINDOWS
#ifdef LOGOG_UNICODE
			nActualSize = vswprintf( pszFormatted, nSizeInWords, cFormatString, argsCopy );
#else // LOGOG_UNICODE
			nActualSize = vsnprintf( pszFormatted, nAttemptedSize, cFormatString, argsCopy );
#endif // LOGOG_UNICODE
#endif // LOGOG_FLAVOR_WINDOWS

			va_end( argsCopy );

			/** Convert the number of LOGOG_CHARs actually formatted into bytes.  This
			 ** does NOT include the trailing NULL.
			 **/
			if ( nActualSize != -1 )
				nActualSize *= sizeof( LOGOG_CHAR );

			/** When we're doing the compare, we have to keep in mind that the nActualSize
			 ** does not include a null.  We need to verify that the nAttemptedSize can hold all
			 ** of nActualSize PLUS the size of one null on this platform.  A LOGOG_CHAR could
			 ** be 1, 2, or 4 bytes long.  So nAttemptedSize must be greater or equal to nActualSize
			 ** less the size of one (null) LOGOG_CHAR in bytes.  Also, the last
			 ** allocation may have failed altogether.
			 ** 
			 **/
			if (( nAttemptedSize >= (nActualSize - (int)sizeof(LOGOG_CHAR))) && ( nActualSize != -1))
				break;

			/** That attempted allocation failed */
			Deallocate( pszFormatted );

			/** If nActualSize has a positive value, it includes the number of bytes needed to hold
			 ** the formatted string; we'll add a LOGOG_CHAR size to the end for the next
			 ** allocation.  If nActualSize has no meaningful value, we'll double the previous
			 ** size and try again.
			 **/
			if (nActualSize > 0)
			{
				nAttemptedSize = nActualSize + sizeof( LOGOG_CHAR );
			}
			else
			{
				nAttemptedSize *= 2;
			}

		}

		m_bIsConst = false;
		assign( pszFormatted );
		/* We just allocated this string, which means it needs to be deallocated
		 * at shutdown time.  The previous function may have changed the const
		 * setting for this string, which means we may need to change it back here... 
		 * */
		m_bIsConst = false;
	}

	void String::Initialize()
	{
		m_pBuffer = NULL;
		m_pOffset = NULL;
		m_pEndOfBuffer = NULL;
		m_pKMP = NULL;
		m_bIsConst = false;
	}

	void String::preKmp( size_t m )
	{
		ScopedLock sl( GetStringSearchMutex() );

		size_t i, j;

		if ( m_pBuffer == NULL )
			return;

		if ( m_pKMP == NULL )
		{
			m_pKMP = (size_t *)Allocate( sizeof( size_t ) * ( m + 1) );
		}

		i = 0;
		j = *m_pKMP = (size_t)-1;

		while (i < m)
		{
			while (j > (size_t)-1 && m_pBuffer[i] != m_pBuffer[j])
				j = m_pKMP[j];
			i++;
			j++;
			if (m_pBuffer[i] == m_pBuffer[j])
				m_pKMP[i] = m_pKMP[j];
			else
				m_pKMP[i] = j;
		}
	}

	size_t String::KMP( const LOGOG_CHAR *y, size_t n )
	{
		size_t i, j;

		size_t m = size() - 1; // ignore NULL char

		/* Preprocessing */
		if ( m_pKMP == NULL )
			preKmp( m );

		/* Searching */
		i = j = 0;
		while (j < n)
		{
			while (i > (size_t)-1 && m_pBuffer[i] != y[j])
				i = m_pKMP[i];
			i++;
			j++;
			if (i >= m)
			{
				return (j - i);
				// We would do this if we cared about multiple substrings
				// i = m_pKMP[i];
			}
		}

		return npos;
	}
}

