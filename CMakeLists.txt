# starter code
cmake_minimum_required(VERSION 3.0)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)

project(Mechanics)

add_library(celestial_mechanics celestial_mechanics.cc)

add_executable(driver driver.cc)
target_link_libraries(driver celestial_mechanics)

# set default build type to Debug
if ( NOT CMAKE_BUILD_TYPE)
  set( CMAKE_BUILD_TYPE Debug CACHE STRING
       "Choose the type of build: Debug Release RelWithDebInfo"
       FORCE)
endif()

add_custom_target(run COMMAND driver)


