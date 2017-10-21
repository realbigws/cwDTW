
aux_source_directory(. MER_SRCS)

MESSAGE(STATUS "building 5mer")

add_library (5mer STATIC ${MER_SRCS})

target_link_libraries(5mer)

set_target_properties(5mer PROPERTIES OUTPUT_NAME "5mer")
set_target_properties(5mer PROPERTIES CLEAN_DIRECT_OUTPUT 1)
