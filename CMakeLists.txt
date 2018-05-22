cmake_minimum_required(VERSION 3.0.0)
project(ice VERSION 0.1.0)

include(CTest)
enable_testing()
aux_source_directory(./src DIRSRCS)
include_directories(./include)

add_subdirectory(./src/alg)
add_executable(ice ${DIRSRCS})
target_link_libraries(ice alg)

install (TARGETS ice DESTINATION ./bin)


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
