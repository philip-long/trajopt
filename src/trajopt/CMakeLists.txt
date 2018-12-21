include_directories(
  SYSTEM
  ${BULLET_INCLUDE_DIRS}  
  ${Boost_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
  ${OpenRAVE_INCLUDE_DIRS}
  ${JSON_INCLUDE_DIRS}
)




if (NOT APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
endif()


#include(CheckCXXCompilerFlag)
#CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
#CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
#if(COMPILER_SUPPORTS_CXX11)
#    set(CMAKE_CXX_FLAGS "-std=c++11")
#elseif(COMPILER_SUPPORTS_CXX0X)
#    set(CMAKE_CXX_FLAGS "-std=c++0x")
#else()
#    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler. Suggested solution: update the pkg build-essential ")
#endif()


add_library(trajopt
	traj_plotter.cpp
  trajectory_costs.cpp
	kinematic_terms.cpp 
	collision_terms.cpp
	json_marshal.cpp
	problem_description.cpp
	bullet_collision_checker.cpp 
	configuration_space.cpp 
	utils.cpp
	rave_utils.cpp
	collision_checker.cpp
	plot_callback.cpp
	bullet_unity.cpp
)
target_link_libraries(trajopt ${OpenRAVE_BOTH_LIBRARIES} ${Boost_SYSTEM_LIBRARY} sco utils json osgviewer)

add_subdirectory(test)

include_directories(${PYTHON_NUMPY_INCLUDE_DIR})
boost_python_module(ctrajoptpy trajoptpy.cpp)
target_link_libraries(ctrajoptpy trajopt ${Boost_SYSTEM_LIBRARY} )

#install(TARGETS trajopt 
#   ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
#)
