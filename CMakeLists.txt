cmake_minimum_required (VERSION 2.6)
project (DOS_CALC)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")

set(BOOST_ROOT "~/BOOST_1_63")
set(BOOST_INCLUDE_DIRS "${BOOST_ROOT}/include")
set(BOOST_LIBRARY_DIRS "${BOOST_ROOT}/lib")

find_package(MPI REQUIRED)
find_package(Boost COMPONENTS program_options mpi serialization REQUIRED)

include_directories("~/BOOST_1_63/include" "/usr/include/openmpi-x86_64")

add_executable(main main.cpp)
target_link_libraries(main ${Boost_LIBRARIES} ${MPI_LIBRARIES})

add_executable(stitch stitch.cpp)
target_link_libraries(stitch ${Boost_LIBRARIES})

add_executable(integrate integrate.cpp)
target_link_libraries(integrate ${Boost_LIBRARIES})
