/**
 * \file macro.hpp Macros for instantiation of a message.
 */

#ifndef __LOGOG_MACRO_HPP__
#define __LOGOG_MACRO_HPP__

namespace logog
{
#ifdef LOGOG_USE_PREFIX
#define LOGOG_PREFIX LOGOG_
#endif // LOGOG_USE_PREFIX

#ifndef LOGOG_GROUP
/** This is the current group for created messages.  Set this to NULL if you
  * want messages to not be part of any specific group.
  */
#define LOGOG_GROUP NULL
#endif

#ifndef LOGOG_CATEGORY
/** This is the current category for created messages.  Set this to NULL if you
  * want messages to not be part of any specific group.
  */
#define LOGOG_CATEGORY NULL
#endif

/** When you have a macro replacement, the preprocessor will only expand the macros recursively
 * if neither the stringizing operator # nor the token-pasting operator ## are applied to it.
 * So, you have to use some extra layers of indirection, you can use the token-pasting operator
 * with a recursively expanded argument.
 */
#define TOKENPASTE2(x, y) x ## y
/** \sa TOKENPASTE2(x, y) */
#define TOKENPASTE(x, y) TOKENPASTE2(x, y)

/** This macro is used when a message is instantiated without any varargs
  * provided by the user.  It locks a global mutex, creates the message,
  * locks it, transmits it, and releases all locks.
  */
#define LOGOG_LEVEL_GROUP_CATEGORY_MESSAGE_NO_VA( level, group, cat, msg ) \
{ \
	Mutex *___pMCM = &GetMessageCreationMutex(); \
	___pMCM->MutexLock(); \
	static logog::Message *TOKENPASTE(_logog_,__LINE__) = new logog::Message( level, \
		LOGOG_CONST_STRING( __FILE__ ), \
		__LINE__ , \
		LOGOG_CONST_STRING( group ), \
		LOGOG_CONST_STRING( cat ), \
		msg; \
	___pMCM->MutexUnlock(); \
	TOKENPASTE(_logog_,__LINE__)->m_Transmitting.MutexLock(); \
	TOKENPASTE(_logog_,__LINE__)->Transmit(); \
	TOKENPASTE(_logog_,__LINE__)->m_Transmitting.MutexUnlock(); \
}

/** This macro is used when a message is instantiated with varargs provided
  * by the user.  It locks a global mutex, creates the message, locks it,
  * formats the message string inside the message, transmits it, 
  * and releases all locks.
  * When logog is shut down, it may be started back up again later.  Therefore,
  * logog needs a way to flag all static Message pointers that they need
  * to be recreated.  We manually simulate a static Message pointer by 
  * implementing it via a static bool.  The bool is turned on the first time
  * this code is run.
  * NOTE!  A subtle race condition exists in the following code, that will ONLY occur
  * if logog is shut down at the same moment that a log message is processed from
  * another thread than the one calling the shutdown.  The Message object could
  * theoretically be destroyed from another thread just before it's locked
  * in a thread that calls the Format() and Transmit() calls on it.  I'm not
  * sure if this is really a bug -- technically, this race condition will
  * occur only if you are calling log messages right on top of the SHUTDOWN
  * call from the main thread.
  */
#define LOGOG_LEVEL_GROUP_CATEGORY_MESSAGE( level, group, cat, formatstring, ... ) \
{ \
	::logog::Mutex *___pMCM = &::logog::GetMessageCreationMutex(); \
	___pMCM->MutexLock(); \
	static bool TOKENPASTE(_logog_static_bool_,__LINE__) = false; \
	static logog::Message * TOKENPASTE(_logog_,__LINE__); \
	if ( TOKENPASTE(_logog_static_bool_,__LINE__) == false ) \
	{ \
		TOKENPASTE(_logog_,__LINE__) = \
			new logog::Message( level, \
				LOGOG_CONST_STRING( __FILE__ ), \
				__LINE__ , \
				LOGOG_CONST_STRING( group ), \
				LOGOG_CONST_STRING( cat ), \
				LOGOG_CONST_STRING( "" ), \
				0.0f, \
				& (TOKENPASTE(_logog_static_bool_,__LINE__)) ); \
	} \
	___pMCM->MutexUnlock(); \
	/* A race condition could theoretically occur here if you are shutting down at the same instant as sending log messages. */ \
	TOKENPASTE(_logog_,__LINE__)->m_Transmitting.MutexLock(); \
	TOKENPASTE(_logog_,__LINE__)->Format( formatstring, ##__VA_ARGS__ ); \
	TOKENPASTE(_logog_,__LINE__)->Transmit(); \
	TOKENPASTE(_logog_,__LINE__)->m_Transmitting.MutexUnlock(); \
}

/** Calls LOGOG_LEVEL_GROUP_CATEGORY_MESSAGE with the current LOGOG_GROUP and
  * LOGOG_CATEGORY setting.
  */
#define LOGOG_LEVEL_MESSAGE( level, formatstring, ... ) \
	LOGOG_LEVEL_GROUP_CATEGORY_MESSAGE( level, LOGOG_GROUP, LOGOG_CATEGORY, formatstring, ##__VA_ARGS__ )

/** Calls LOGOG_LEVEL_MESSAGE with the current LOGOG_LEVEL setting. */
#define LOGOG_MESSAGE( formatstring, ... ) \
	LOGOG_LEVEL_MESSAGE( LOGOG_LEVEL, formatstring, ##__VA_ARGS__ )


#if LOGOG_LEVEL >= LOGOG_LEVEL_DEBUG
/** Logs a message at the DEBUG reporting level. */
#define LOGOG_DEBUG( formatstring, ... ) \
	LOGOG_LEVEL_MESSAGE( LOGOG_LEVEL_DEBUG, formatstring, ##__VA_ARGS__ )
#else
#define LOGOG_DEBUG( formatstring, ... ) {};
#endif

#if LOGOG_LEVEL >= LOGOG_LEVEL_INFO
/** Logs a message at the INFO reporting level. */
#define LOGOG_INFO( formatstring, ... ) \
	LOGOG_LEVEL_MESSAGE( LOGOG_LEVEL_INFO, formatstring, ##__VA_ARGS__ )
#else
#define LOGOG_INFO( formatstring, ... ) {};
#endif

#if LOGOG_LEVEL	>= LOGOG_LEVEL_WARN3
/** Logs a message at the WARN3 reporting level. */
#define LOGOG_WARN3( formatstring, ... ) \
	LOGOG_LEVEL_MESSAGE( LOGOG_LEVEL_WARN3, formatstring, ##__VA_ARGS__ )
#else
#define LOGOG_WARN3( formatstring, ... ) {};
#endif

#if LOGOG_LEVEL	>= LOGOG_LEVEL_WARN2
/** Logs a message at the WARN2 reporting level. */
#define LOGOG_WARN2( formatstring, ... ) \
	LOGOG_LEVEL_MESSAGE( LOGOG_LEVEL_WARN2, formatstring, ##__VA_ARGS__ )
#else
#define LOGOG_WARN2( formatstring, ... ) {};
#endif

#if LOGOG_LEVEL	>= LOGOG_LEVEL_WARN1
/** Logs a message at the WARN1 reporting level. */
#define LOGOG_WARN1( formatstring, ... ) \
	LOGOG_LEVEL_MESSAGE( LOGOG_LEVEL_WARN1, formatstring, ##__VA_ARGS__ )
#else
#define LOGOG_WARN1( formatstring, ... ) {};
#endif

#if LOGOG_LEVEL	>= LOGOG_LEVEL_WARN
/** Logs a message at the WARN reporting level. */
#define LOGOG_WARN( formatstring, ... ) \
	LOGOG_LEVEL_MESSAGE( LOGOG_LEVEL_WARN, formatstring, ##__VA_ARGS__ )
#else
#define LOGOG_WARN( formatstring, ... ) {};
#endif

#if LOGOG_LEVEL	>= LOGOG_LEVEL_ERROR
/** Logs a message at the ERROR reporting level. */
#define LOGOG_ERROR( formatstring, ... ) \
	LOGOG_LEVEL_MESSAGE( LOGOG_LEVEL_ERROR, formatstring, ##__VA_ARGS__ )
#else
#define LOGOG_ERROR( formatstring, ... ) {};
#endif

#if LOGOG_LEVEL	>= LOGOG_LEVEL_CRITICAL
/** Logs a message at the CRITICAL reporting level. */
#define LOGOG_CRITICAL( formatstring, ... ) \
	LOGOG_LEVEL_MESSAGE( LOGOG_LEVEL_CRITICAL, formatstring, ##__VA_ARGS__ )
#else
#define LOGOG_CRITICAL( formatstring, ... ) {};
#endif

#if LOGOG_LEVEL	>= LOGOG_LEVEL_ALERT
/** Logs a message at the ALERT reporting level. */
#define LOGOG_ALERT( formatstring, ... ) \
	LOGOG_LEVEL_MESSAGE( LOGOG_LEVEL_ALERT, formatstring, ##__VA_ARGS__ )
#else
#define LOGOG_ALERT( formatstring, ... ) {};
#endif

#if LOGOG_LEVEL	>= LOGOG_LEVEL_EMERGENCY
/** Logs a message at the EMERGENCY reporting level. */
#define LOGOG_EMERGENCY( formatstring, ... ) \
	LOGOG_LEVEL_MESSAGE( LOGOG_LEVEL_EMERGENCY, formatstring, ##__VA_ARGS__ )
#else
#define LOGOG_EMERGENCY( formatstring, ... ) {};
#endif

#define LOGOG_SET_LEVEL( level ) \
	::logog::SetDefaultLevel( level );

/** Define this compilation flag if your compilation environment conflicts with
  * any of the shorthand logging macros in macro.hpp.
  */
#ifndef LOGOG_USE_PREFIX
/* If you get compilation errors in this section, then define the flag LOGOG_USE_PREFIX during compilation, and these
 * shorthand logging macros won't exist -- you'll need to use the LOGOG_* equivalents above.
 */
/* We can't use DEBUG in Win32 unfortunately, so we use DBUG for shorthand here. */
//! [Shorthand]
/** \sa LOGOG_DEBUG */
#define DBUG(...) LOGOG_DEBUG( __VA_ARGS__ )
/** \sa LOGOG_INFO */
#define INFO(...) LOGOG_INFO( __VA_ARGS__ )
/** \sa LOGOG_WARN3 */
#define WARN3(...) LOGOG_WARN3( __VA_ARGS__ )
/** \sa LOGOG_WARN2 */
#define WARN2(...) LOGOG_WARN2( __VA_ARGS__ )
/** \sa LOGOG_WARN1 */
#define WARN1(...) LOGOG_WARN1( __VA_ARGS__ )
/** \sa LOGOG_WARN */
#define WARN(...) LOGOG_WARN( __VA_ARGS__ )
/** \sa LOGOG_ERROR */
#define ERR(...) LOGOG_ERROR( __VA_ARGS__ )
/** \sa LOGOG_ALERT */
#define ALERT(...) LOGOG_ALERT( __VA_ARGS__ )
/** \sa LOGOG_CRITICAL */
#define CRITICAL(...) LOGOG_CRITICAL( __VA_ARGS__ )
/** \sa LOGOG_EMERGENCY */
#define EMERGENCY(...) LOGOG_EMERGENCY( __VA_ARGS__ )
//! [Shorthand]
#endif

/** Call this function to initialize logog and prepare for logging. 
  * \sa logog::Initialize()
  */
#define LOGOG_INITIALIZE(...)  logog::Initialize( __VA_ARGS__ );

/** Call this function to shut down logog and release all memory allocated.
  * \sa logog::Shutdown()
  */

#define LOGOG_SHUTDOWN()   logog::Shutdown();

} // namespace logog

#endif // __LOGOG_MACRO_HPP_
