
aux_source_directory(. MER_SRCS)

MESSAGE(STATUS "building kmer")

add_library (kmer STATIC ${MER_SRCS})

target_link_libraries(kmer)

set_target_properties(kmer PROPERTIES OUTPUT_NAME "kmer")
set_target_properties(kmer PROPERTIES CLEAN_DIRECT_OUTPUT 1)
