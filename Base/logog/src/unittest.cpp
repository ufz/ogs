
/*
 * \file unittest.cpp
 */

#define LOGOG_UNIT_TESTING 1

#include "logog.hpp"

namespace logog
{

TestRegistryType &LogogTestRegistry()
{
    static TestRegistryType *pRegistry = new TestRegistryType();
    return *pRegistry;
}

TestSignup::TestSignup( UnitTest *pTest )
{
    m_pTest = pTest;
    LogogTestRegistry().push_back( pTest );
}

UnitTest::UnitTest( const TestNameType &sTestName )
{
    m_sTestName = sTestName;
    m_pTestSignup = new TestSignup( this );
}

UnitTest::~UnitTest()
{
	FreeInternals();
}

void UnitTest::FreeInternals()
{
	if ( m_pTestSignup )
		delete m_pTestSignup;

	m_pTestSignup = NULL;
}

/** Returns the name of this UnitTest provided at construction time. */
TestNameType &UnitTest::GetName()
{
	return m_sTestName;
}

/** Executes all currently registered tests and prints a report of success or failure. */
int RunAllTests()
{
    using namespace std;

    int nTests = 0, nTestsSucceeded = 0;
    int nTestResult;
    int nFailures = 0;

#ifdef LOGOG_UNICODE
    wostream *pOut;
#else // LOGOG_UNICODE
	ostream *pOut;	
#endif // LOGOG_UNICODE

	pOut = &(LOGOG_COUT);

    nTests = (int) LogogTestRegistry().size();

    if ( nTests == 0 )
    {
		*pOut << _LG("No tests currently defined.") << endl;
        return 1;
    }

    for ( TestRegistryType::iterator it = LogogTestRegistry().begin();
            it != LogogTestRegistry().end();
            ++it )
    {
        (*pOut) << _LG("Test ") << (*it)->GetName() << _LG(" running... ") << endl;
        nTestResult = (*it)->RunTest();

        (*pOut) << _LG("Test ") << (*it)->GetName();

        if ( nTestResult == 0 )
        {
            *pOut << _LG(" successful.") << endl;
            nTestsSucceeded++;
        }
        else
        {
            *pOut << _LG(" failed!") << endl;
            nFailures++;
        }

        /* Validate that no allocations are currently outstanding.  Make sure to handle the case
         * where leak detection is disabled */
        int nMemoryTestResult = ReportMemoryAllocations();

        if ( nMemoryTestResult != -1 )
        {
            (*pOut) << _LG("Test ") << (*it)->GetName() << _LG(" has ") << nMemoryTestResult <<
				_LG(" memory allocations outstanding at end of test.") << endl;
            nFailures += nMemoryTestResult;
        }
    }

    *pOut << _LG("Testing complete, ")
          << nTests << _LG(" total tests, ")
          << nTestsSucceeded << _LG(" tests succeeded, ")
          << ( nTests - nTestsSucceeded ) << _LG(" failed")
          << endl;

    return nFailures;
}

/** Should remove all memory allocated during unit testing. */
void ShutdownTests()
{
	TestRegistryType::iterator it;

	for ( it = LogogTestRegistry().begin(); it != LogogTestRegistry().end(); it++ )
		(*it)->FreeInternals();

    delete &LogogTestRegistry();
}

}

