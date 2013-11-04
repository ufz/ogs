/**
 * \file platform.hpp Platform specific detection and defines, including STL overrides
 */

#ifndef __LOGOG_PLATFORM_HPP__
#define __LOGOG_PLATFORM_HPP__

#ifdef DOXYGEN
/** Set this compilation flag to 1 to indicate that this is a Windows-like platform, e.g. Win32, Win64, Xbox 360, etc. */
#define LOGOG_FLAVOR_WINDOWS 1
/**  Set this compilation flag to 1 to indicate that this is a Posix-like platform, e.g. Linux, PS3, etc. */
#define LOGOG_FLAVOR_POSIX 1
#endif

#ifndef LOGOG_FLAVOR_WINDOWS
#ifndef LOGOG_FLAVOR_POSIX

/* The user hasn't told us which flavor we're running on, so we must make a guess as to which platform is valid. */

/* If this is MSVC, then it's Windows like */
#ifdef _MSC_VER
#define LOGOG_FLAVOR_WINDOWS 1
#endif // _MSC_VER

/* gcc probably means Posix */
#ifdef __GNUC__
#define LOGOG_FLAVOR_POSIX 1
#endif

#ifdef __CYGWIN__
#ifdef LOGOG_UNICODE
#error LOGOG_UNICODE not supported on Cygwin platform because Cygwin does not support vsnwprintf
#endif
/* Cygwin lacks vsnprintf support in headers but it's in the libraries.  First we
 * need to define the constants size_t and va_list, then define vsnprintf */
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
extern int vsnprintf(char *str, size_t size, const char *format, va_list ap);
#define LOGOG_USE_TR1 1
#endif // __CYGWIN__

#ifdef __linux__
#define LOGOG_USE_TR1 1
#endif

#ifdef __APPLE__
#ifndef __clang__
#define LOGOG_USE_TR1 1
#else
#if !(__has_include(<unordered_map>))
#define LOGOG_USE_TR1 1
#endif // !(__has_include(<unordered_map>))
#endif // __clang
#endif // __APPLE__

/* Detect IBM's XL C++ */
#ifdef __IBMCPP__
// Enable use of TR1 unorderd_map etc.
#define __IBMCPP_TR1__ 1
#ifdef __PPC__
#define LOGOG_FLAVOR_POSIX 1
#endif // __PPC__
#endif // __IBMCPP__

/* If we've recognized it already, it's a relatively modern compiler */
#if defined( LOGOG_FLAVOR_WINDOWS ) || defined( LOGOG_FLAVOR_POSIX )
#define LOGOG_HAS_UNORDERED_MAP 1
#endif // defined(...)

/* PS3 */
#ifdef SN_TARGET_PS3
#include "proprietary/ps3.hpp"
#endif

#endif // LOGOG_FLAVOR_POSIX
#endif // LOGOG_FLAVOR_WINDOWS

#ifdef LOGOG_FLAVOR_WINDOWS
/* Detect Xbox 360 */
#if _XBOX_VER >= 200
#include "proprietary/xbox360.hpp"
#else
/* Windows has been detected. */

/** Microsoft's CRT library has its own leak detection mechanism.  If you don't trust logog's
 ** LOGOG_LEAK_DETECTION, you can enable LOGOG_LEAK_DETECTION_WINDOWS to enable
 ** Microsoft's version.  This really doesn't belong in platform.hpp -- feel
 ** free to refactor this into a Windows-specific header.
 **/
#ifdef LOGOG_LEAK_DETECTION_WINDOWS
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

#ifdef _DEBUG
#ifndef DBG_NEW
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
#define new DBG_NEW
#endif
#endif  // _DEBUG
#endif // LOGOG_LEAK_DETECTION_WINDOWS

#include "windows.h"
#endif // _XBOX_VER
#endif // LOGOG_FLAVOR_WINDOWS

#ifndef LOGOG_FLAVOR_WINDOWS
#ifndef LOGOG_FLAVOR_POSIX
#error Platform flavor not detected.  Please see platform.hpp to configure this platform.
#endif
#endif

// For the FILE type.
#include <cstdio>
// For POSIX file access.
#include <cstdlib>
// For Ubuntu 11 and ptrdiff_t
#include <cstddef>

#include <ctime>

#if defined(LOGOG_LEAK_DETECTION) || defined(LOGOG_INTERNAL_DEBUGGING)
#include <cstdio>
#include <iostream>
#endif // LOGOG_LEAK_DETECTION || LOGOG_INTERNAL_DEBUGGING

#ifdef LOGOG_FLAVOR_POSIX
#include <sys/time.h>
// for uint64_t
#include <inttypes.h>
#include <wchar.h>
// for strstr support
#include <cstring>
#endif

// For Unicode support
#include <typeinfo>

/* ----------------------------------------------------------- */
/* Here's the stuff your compiler may have a problem with...   */

//! [STLTypes]

/** The definition for a hash type in logog.  You can replace it here with your own set compatible class if needed. */
#define LOGOG_HASH			std::tr1::hash
/** The definition for an unordered map type in logog.  You can replace it here with your own unordered_map compatible class if needed. */

#ifdef LOGOG_HAS_UNORDERED_MAP
#ifdef LOGOG_USE_TR1
#include <tr1/unordered_map>
#else
#include <unordered_map>
#endif // LOGOG_USE_TR1
#ifdef _GLIBCXX_UNORDERED_MAP
#define LOGOG_UNORDERED_MAP	std::unordered_map
#else // _GLIBCXX_UNORDERED_MAP
#define LOGOG_UNORDERED_MAP	std::tr1::unordered_map
#endif  // _GLIBCXX_UNORDERED_MAP
#else // LOGOG_HAS_UNORDERED_MAP
#include <hash_map>
#define LOGOG_UNORDERED_MAP std::hash_map
#endif // LOGOG_HAS_UNORDERED_MAP

/** An internal consistency error has been detected in logog.  */
#define LOGOG_INTERNAL_FAILURE	abort();

/* ----------------------------------------------- */
/* Here's the stuff that's pretty standard in STL by now. */

/** The definition for a pair type in logog.  You can replace it here with your own pair compatible class if needed. */
#define LOGOG_PAIR			std::pair
#include <list>
/** The definition for a list type in logog.  You can replace it here with your own list compatible class if needed. */
#define LOGOG_LIST			std::list
#include <vector>
/** The definition for a vector type in logog.  You can replace it here with your own vector compatible class if needed. */
#define LOGOG_VECTOR		std::vector
#include <set>
/** The definition for a set type in logog.  You can replace it here with your own set compatible class if needed. */
#define LOGOG_SET			std::set
/** The definition for STL "equal to" in logog.  You can replace it here with your own set compatible class if needed. */
#define LOGOG_EQUAL_TO		std::equal_to

//! [STLTypes]

/** The default port number that logog uses to communicate via TCP/UDP */
#define LOGOG_DEFAULT_PORT	9987

/** We use va_list for formatting messages. */
#include <cstdarg>



#endif // __LOGOG_PLATFORM_HPP
