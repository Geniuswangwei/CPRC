cmake_minimum_required(VERSION 3.10)
project(honeycomb_solver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add compile options
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
endif()

# Find OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Add demo_1 executable
add_executable(demo_1 demo_1.cpp)

# Set include directories for demo_1
target_include_directories(demo_1 PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link OpenMP for demo_1
target_link_libraries(demo_1 PUBLIC 
    OpenMP::OpenMP_CXX
)

# Install
install(TARGETS demo_1
    RUNTIME DESTINATION bin
)