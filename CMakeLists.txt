project(statistics-cgp)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -std=c++11)

set(CMAKE_MODULE_PATH "${statistics-cgp_SOURCE_DIR}/config")

find_package(Boost COMPONENTS program_options REQUIRED)
find_package(HDF5 COMPONENTS C HL REQUIRED)
find_package(CGP REQUIRED)
find_package(Chatty REQUIRED)
find_package(VecVec REQUIRED)
find_package(Vigra REQUIRED)

include_directories(${CGP_INCLUDE_DIR})

add_executable(cgp_statistics
    compressors.cxx
    supervoxels.cxx
    cgp_statistics.cxx
)
target_link_libraries(cgp_statistics
    ${HDF5_LIBRARIES}
    ${VIGRA_IMPEX_LIBRARY}
    ${CHATTY_LIBRARY}
    ${VECVEC_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
)

get_property(location TARGET cgp_statistics PROPERTY LOCATION)
add_custom_command(TARGET cgp_statistics
    POST_BUILD        
    COMMAND ${CMAKE_COMMAND} -E copy_if_different  
        ${location}
        ${CMAKE_CURRENT_SOURCE_DIR}
)

