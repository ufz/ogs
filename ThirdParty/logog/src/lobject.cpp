/* 
 * \file object.cpp
 */

#include "logog.hpp"
#ifdef LOGOG_LEAK_DETECTION
#include <iostream>
#endif // LOGOG_LEAK_DETECTION

namespace logog {

#ifdef LOGOG_LEAK_DETECTION
AllocationsType s_Allocations;
#endif

	Object::Object() {}

	Object::~Object()
	{

	}

	void *Object::operator new( size_t nSize )
	{
		return Allocate( nSize );
	}

	void *Object::operator new[]( size_t nSize )
	{
		return Allocate( nSize );
	}

#ifdef _DEBUG
#ifdef LOGOG_FLAVOR_WINDOWS

	void *Object::operator new(size_t nSize, LPCSTR lpszFileName, int nLine)
	{
		/* avoid unref'd parameter warnings */
		lpszFileName;
		nLine;
		return Allocate( nSize );
	}

	void *Object::operator new[](size_t nSize, LPCSTR lpszFileName, int nLine)
	{
		/* avoid unref'd parameter warnings */
		lpszFileName;
		nLine;
		return Allocate( nSize );
	}

	void Object::operator delete(void *ptr, LPCSTR lpszFileName, int nLine)
	{
		lpszFileName;
		nLine;
		Deallocate( ptr );
	}

	void Object::operator delete[](void *ptr, LPCSTR lpszFileName, int nLine)
	{
		lpszFileName;
		nLine;
		Deallocate( ptr );
	}

#endif // LOGOG_FLAVOR_WINDOWS
#endif // _DEBUG

	    /** Deletes an object pointed to by ptr. */
	void Object::operator delete( void *ptr )
    {
        Deallocate( ptr );
    }
    /** Deletes an object array pointed to by ptr. */
	void Object::operator delete[]( void *ptr )
    {
        Deallocate( ptr );
    }

    /** Allocates nSize bytes of memory.  You must call logog::Initialize() before calling this function.
     * \sa Initialize()
     */
	void *Object::Allocate( size_t nSize )
    {
        void *ptr = Static().s_pfMalloc( nSize );
#ifdef LOGOG_REPORT_ALLOCATIONS
        LOGOG_COUT << _LG("Allocated ") << nSize << _LG(" bytes of memory at ") << ptr << endl;
#endif // LOGOG_REPORT_ALLOCATIONS
#ifdef LOGOG_LEAK_DETECTION
        AllocationsType::iterator it;

        LockAllocationsMutex();
        it = s_Allocations.find( ptr );

        if ( it != s_Allocations.end() )
        {
            LOGOG_COUT << _LG("Reallocation detected in memory manager!  We seem to have allocated the same address twice ")
                 << _LG("without freeing it!  Address = ") << ptr << std::endl;
            UnlockAllocationsMutex();
            LOGOG_INTERNAL_FAILURE;
        }

        s_Allocations.insert( LOGOG_PAIR< const PointerType, size_t >( ptr, nSize ) );
        UnlockAllocationsMutex();
#endif // LOGOG_LEAK_DETECTION
        return ptr;
    }

    /** Deallocate a pointer previously acquired by Allocate(). */
	void Object::Deallocate( void *ptr )
    {
#ifdef LOGOG_LEAK_DETECTION
        LockAllocationsMutex();
        AllocationsType::iterator it;

        it = s_Allocations.find( ptr );

        if ( it == s_Allocations.end() )
        {
            LOGOG_COUT << _LG("Freeing memory not previously allocated!  Address = ") << ptr << std::endl;
            UnlockAllocationsMutex();
            LOGOG_INTERNAL_FAILURE;
        }

#ifdef LOGOG_REPORT_ALLOCATIONS
        LOGOG_COUT << _LG("Freeing ") << it->second << _LG(" bytes of memory at ") << it->first << endl;
#endif // LOGOG_REPORT_ALLOCATIONS
        s_Allocations.erase( ptr );
        UnlockAllocationsMutex();
#endif // LOGOG_LEAK_DETECTION
        Static().s_pfFree( ptr );
    }

	int MemoryAllocations()
	{
#ifdef LOGOG_LEAK_DETECTION
		LockAllocationsMutex();
		size_t nSize = s_Allocations.size();

		if ( nSize != 0 )
			LOGOG_COUT << _LG("Total active allocations: ") << nSize << std::endl;

		UnlockAllocationsMutex();
		return (int)nSize;
#else // LOGOG_LEAK_DETECTION
		return -1;
#endif // LOGOG_LEAK_DETECTION
	}

	int ReportMemoryAllocations()
	{
#ifdef LOGOG_LEAK_DETECTION
		LockAllocationsMutex();

		if ( s_Allocations.size() == 0 )
		{
			LOGOG_COUT << _LG("No memory allocations outstanding.") << std::endl;
		}
		else
		{
			for ( AllocationsType::iterator it = s_Allocations.begin();
				it != s_Allocations.end();
				it++ )
			{
				LOGOG_COUT << _LG("Memory allocated at ") << it->first << 
					_LG(" with size ") << it->second << 
					_LG(" bytes ") << std::endl;
			}
		}

		UnlockAllocationsMutex();
#endif
		return MemoryAllocations();
	}

}

