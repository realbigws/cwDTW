
aux_source_directory(. PROC_SRCS)

MESSAGE(STATUS "building proc")

include_directories(${PROJECT_SOURCE_DIR})

add_library (proc STATIC ${PROC_SRCS})

target_link_libraries(proc)

set_target_properties(proc PROPERTIES OUTPUT_NAME "proc")
set_target_properties(proc PROPERTIES CLEAN_DIRECT_OUTPUT 1)
