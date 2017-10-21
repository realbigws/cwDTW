
MESSAGE(STATUS "building main")

add_executable(wletdtw main.cpp)
target_link_libraries(wletdtw proc 5mer m wavelib)
