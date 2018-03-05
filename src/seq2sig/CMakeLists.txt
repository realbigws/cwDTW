
MESSAGE(STATUS "building seq2sig")

add_executable(seq2sig main.cpp)

target_link_libraries(seq2sig kmer proc m)
