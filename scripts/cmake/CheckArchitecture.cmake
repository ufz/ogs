# Check architecture
if(MSVC AND NOT HAVE_64_BIT AND NOT OGS_32_BIT)
    message(FATAL_ERROR "Building OGS on Windows with 32-bit is not supported! \
Either use the correct generator, e.g. 'Visual Studio 14 2015 Win64' or define \
'-DOGS_32_BIT=ON' if you know what you are doing.")
endif()
