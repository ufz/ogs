APPEND_SOURCE_FILES(SOURCES_CTIVE_REL)

message(WARNING "sources ${SOURCES_CTIVE_REL}")

add_library(ConstitutiveRelationsLib ${SOURCES_CTIVE_REL})

target_link_libraries(ConstitutiveRelationsLib ${VTK_LIBRARIES})

ADD_VTK_DEPENDENCY(ConstitutiveRelationsLib)
