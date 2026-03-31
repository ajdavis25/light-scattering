cmake_minimum_required(VERSION 3.10)
project(MonteCarloCPP)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_USE_RESPONSE_FILE_FOR_INCLUDES OFF)

enable_testing()
find_package(OpenMP)

add_library(rt_core
    src/Atmosphere.cpp
    src/WavelengthHandling.cpp
    src/PhaseFunctions.cpp
    src/Polarization.cpp
    src/SurfaceReflection.cpp
    src/MonteCarloDriver.cpp
    src/ValidationQA.cpp
)

target_include_directories(rt_core PUBLIC src)
if(OpenMP_CXX_FOUND)
    target_link_libraries(rt_core PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(MonteCarloCPP src/main.cpp)
target_link_libraries(MonteCarloCPP PRIVATE rt_core)

add_executable(ValidationRunner src/validation_main.cpp)
target_link_libraries(ValidationRunner PRIVATE rt_core)

add_executable(MeasurementCaseRunner src/measurement_case_main.cpp)
target_link_libraries(MeasurementCaseRunner PRIVATE rt_core)

add_executable(test_phasefunctions tests/test_phasefunctions.cpp)
target_link_libraries(test_phasefunctions PRIVATE rt_core)

add_executable(test_surface tests/test_surface.cpp)
target_link_libraries(test_surface PRIVATE rt_core)

add_executable(test_polarization tests/test_polarization.cpp)
target_link_libraries(test_polarization PRIVATE rt_core)

add_executable(test_atmosphere tests/test_atmosphere.cpp)
target_link_libraries(test_atmosphere PRIVATE rt_core)

add_executable(test_wavelength tests/test_wavelength.cpp)
target_link_libraries(test_wavelength PRIVATE rt_core)

add_test(NAME phasefunctions COMMAND test_phasefunctions)
add_test(NAME surface COMMAND test_surface)
add_test(NAME polarization COMMAND test_polarization)
add_test(NAME atmosphere COMMAND test_atmosphere)
add_test(NAME wavelength COMMAND test_wavelength)
