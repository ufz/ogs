include (CheckTypeSize)

CHECK_TYPE_SIZE(int SIZEOF_INT)
CHECK_TYPE_SIZE(long SIZEOF_LONG)
CHECK_TYPE_SIZE("long long" SIZEOF_LONG_LONG)
CHECK_TYPE_SIZE("void *" SIZEOF_VOID_P)

# Sets Bits variables
# check 64 bit
if( CMAKE_SIZEOF_VOID_P EQUAL 4 )
	set( HAVE_64_BIT 0 )
	set( BITS 32 )
else()
	set( HAVE_64_BIT 1 )
	add_definitions(-DHAVE_64_BIT)
	set( BITS 64)
endif()

# Check endian of the system
include (TestBigEndian)
TEST_BIG_ENDIAN (IS_BIG_ENDIAN)
if (IS_BIG_ENDIAN)
    add_definitions(-DENDIAN_IS_BIG)
   # cannot use BIG_ENDIAN because it's reserved in Linux
else ()
    add_definitions(-DENDIAN_IS_LITTLE)
   # cannot use LITTLE_ENDIAN because it's reserved in Linux
endif () # IS_BIG_ENDIAN

