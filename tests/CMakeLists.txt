cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(Catch2)

# Local
set(SOURCE_FILES
	sample_tests.cpp)
    
add_executable(learnSPH_tests ${SOURCE_FILES})
target_include_directories(learnSPH_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../learnSPH/include)
target_link_libraries(learnSPH_tests learnSPH Catch2::Catch2WithMain)


#add_dependencies(learnSPH_tests learnSPH)
