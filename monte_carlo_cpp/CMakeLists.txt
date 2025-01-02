# CMakeLists.txt
cmake_minimum_required(VERSION 3.0)
project(MonteCarloCPP)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create an executable from these sources
add_executable(MonteCarloCPP
    src/main.cpp
    src/PhaseFunctions.cpp
    src/SurfaceReflection.cpp
    src/Polarization.cpp
    src/Atmosphere.cpp
    src/WavelengthHandling.cpp
    src/AccelerationTechniques.cpp
    # etc...
)
