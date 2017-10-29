
MESSAGE(STATUS "building sig2peak")

add_executable(sig2peak main.cpp)
target_link_libraries(sig2peak proc m wavelib)


