add_library ( swmm5interface
    swmm5_iface.h
    swmm5_iface.c
)

target_link_libraries (swmm5interface swmm5)

include_directories ( swmm5 )
