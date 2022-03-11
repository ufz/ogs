# From https://www.youtube.com/watch?v=8y7UuAG3Z0g (minute 52)
cmake_host_system_information(RESULT _memfree QUERY AVAILABLE_PHYSICAL_MEMORY)
cmake_host_system_information(RESULT _cores QUERY NUMBER_OF_PHYSICAL_CORES)
message(
    STATUS "Number of (logical) cores: ${_cores}, Free memory: ${_memfree} MB"
)

# Sets number of jobs between 1 and number of logical cores depending on the
# available memory.
function(setup_job_pool name mem_per_task)
    math(EXPR res "${_memfree} / ${mem_per_task}")
    if(res LESS 1)
        set(res 1)
    endif()
    if(res GREATER ${_cores})
        set(res ${_cores})
    endif()
    message(STATUS "  Job pool ${name} using ${res} cores.")
    set_property(GLOBAL APPEND PROPERTY JOB_POOLS ${name}=${res})
endfunction()

# Default job pool
setup_job_pool(light_tasks 800) # MB per task
set(CMAKE_JOB_POOL_COMPILE light_tasks)
set(CMAKE_JOB_POOL_LINK light_tasks)

if(APPLE_ARM)
    setup_job_pool(heavy_tasks 2500)
else()
    setup_job_pool(heavy_tasks 6000)
endif()
