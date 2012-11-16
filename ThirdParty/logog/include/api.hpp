/**
 * \file api.hpp Initialization and shutdown functions for logog.
 */

#ifndef __LOGOG_API_HPP__
#define __LOGOG_API_HPP__

namespace logog
{
//! [INIT_PARAMS]
/**
 * To initialize the memory manager with non-default values, allocate a temporary INIT_PARAMS structure, fill it with
 * zeros, and then change the individual entries in the INIT_PARAMS struct before passing as a parameter to Initialize().
 * \sa Initialize
 */
struct INIT_PARAMS
{
    /** A pointer to a function that allocates memory.  By default logog allocates using malloc().  Change this pointer
     ** before passing to Initialize() to use your own custom memory allocator.
     * \sa logog::Initialize()
     */
    void *( *m_pfMalloc )( size_t );

    /** A pointer to a function that frees memory.  By default logog frees using free().  Change this pointer
     * before passing to Initialize() to use your own custom memory allocator.
     * \sa logog::Initialize()
     */
    void ( *m_pfFree )( void * );
};
//! [INIT_PARAMS]

/** Initializes the logog system.  No logog calls or allocations may be made before calling this function; expect
 * crashes if you haven't called this at the top of your program.
 * \param params The address of an INIT_PARAMS structure you have already allocated on the heap, or NULL to
 * use default values.
 * \sa INIT_PARAMS
 */
extern int Initialize( INIT_PARAMS *params = NULL );

/** Shuts down the logog system and frees all memory allocated by logog.  Memory still allocated by the logog system after Shutdown() indicates
 ** a bug.
 **/
extern int Shutdown( );

}

#endif // __LOGOG_API_HPP
