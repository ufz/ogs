/**
 * \file string.hpp Defines the logog string class.
 */
#ifndef __LOGOG_STRING_HPP__
#define __LOGOG_STRING_HPP__

/** \def LOGOG_UNICODE
 ** Define this typedef in order to support wide characters in all logog strings.
 ** \sa _LG
 ** */

#ifdef LOGOG_UNICODE
/** logog has detected Unicode; therefore a LOGOG_CHAR is a wide character. */
typedef wchar_t LOGOG_CHAR;
/** This is a naughty little hack.  If we're using Unicode, then we append an
 * L to your const string.  However, we have several situations in which that string will actually be a NULL or 0
 * value in code, and the string will be rendered as L0 or LNULL.  In that case, we catch it with yet another 
 * macro.  Hacky, but seems to do the trick.  Beware of conflicts with existing code though...!
 */
#define L0 (const LOGOG_CHAR *)'\0'
#define LNULL (const LOGOG_CHAR *)'\0'
#define L__null (const LOGOG_CHAR *)'\0'

/* Two indirect references are necessary in order to expand a string and append the L to it. */
#define LOGOG_CONST_STRING_INDIRECT(x) L ## x
/* This macro will cause a const string to be stored as a Unicode string on a Unicode build. */
#define LOGOG_CONST_STRING(x) LOGOG_CONST_STRING_INDIRECT(x)
#define LOGOG_COUT  std::wcout
#define LOGOG_CERR  std::wcerr

#else // LOGOG_UNICODE

/** logog has not detected Unicode; therefore a LOGOG_CHAR is simply a char. */
typedef char LOGOG_CHAR;
/** This macro will cause a const string to be stored as an ANSI string on an ANSI build, and as a Unicode string
 * on a Unicode build.
 */
#define LOGOG_CONST_STRING(x) (x)
#define LOGOG_COUT  std::cout
#define LOGOG_CERR  std::cerr
#endif // LOGOG_UNICODE

/** If this constant is defined, then you can use the shorthand macro _LG in your code to represent a constant
  * string.
  */
#ifndef LOGOG_USE_PREFIX
/** The _LG() macro is defined only if LOGOG_USE_PREFIX is not defined.  _LG() can be used to describe 
 ** a const string that is compiled to either as Unicode
 * or ANSI, based on the setting of the LOGOG_UNICODE flag.
 * _LG() is not needed if you don't need Unicode support.  If you want your messages to work with both Unicode 
 * as well as ANSI builds of logog, preface them like this: _LG("This const string works on Unicode as well as ANSI.")
 */
#define _LG( x ) LOGOG_CONST_STRING( x )
#endif 

namespace logog
{
	class String : public Object
	{
	public:

		static const size_t npos = (size_t) -1;

		String();
		virtual ~String();
		virtual void Free();
		static size_t Length( const LOGOG_CHAR *chars );

		String( const String &other );
		String( const LOGOG_CHAR *pstr );
		String & operator =( const String & other);
		String & operator =( const LOGOG_CHAR *pstr );
		size_t size() const;
		virtual void clear();
		virtual size_t reserve( size_t nSize );
		virtual size_t reserve_for_int();
		virtual operator const LOGOG_CHAR *() const;
		virtual size_t assign( const String &other );
		virtual size_t append( const String &other );
		virtual size_t append( const LOGOG_CHAR *other );
		virtual void reverse( LOGOG_CHAR* pStart, LOGOG_CHAR* pEnd);
		virtual size_t assign( const int value );
		virtual size_t append( const LOGOG_CHAR c );
		virtual bool is_valid() const;
		virtual size_t assign( const LOGOG_CHAR *other, const LOGOG_CHAR *pEnd = NULL );

		virtual size_t find( String &other ) const;
		virtual void format( const LOGOG_CHAR *cFormatString, ... );
		virtual void format_va( const LOGOG_CHAR *cFormatString, va_list args );

		virtual const LOGOG_CHAR* c_str() const;

	protected:
		virtual void Initialize();

		/* Code modified from http://www-igm.univ-mlv.fr/~lecroq/string/node8.html#SECTION0080 */
		void preKmp(size_t m);

		size_t KMP( const LOGOG_CHAR *y, size_t n );

#define LOGOG_MAX( a, b ) ( ( a > b ) ? a : b )

		size_t BM(LOGOG_CHAR *y, size_t n);

		LOGOG_CHAR *m_pBuffer;
		LOGOG_CHAR *m_pOffset;
		LOGOG_CHAR *m_pEndOfBuffer;
		size_t *m_pKMP;
		bool m_bIsConst;
	};
}

#define LOGOG_STRING ::logog::String

#endif // __LOGOG_STRING_HPP_
