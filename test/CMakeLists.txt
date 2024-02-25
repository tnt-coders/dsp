find_package(Catch2 REQUIRED)
find_package(FFTW3 REQUIRED)

add_executable(${PROJECT_NAME}_test
    main.cpp
    fourier_transform.cpp
)

target_compile_options(${PROJECT_NAME}_test PRIVATE -D_USE_MATH_DEFINES)

target_link_libraries(${PROJECT_NAME}_test PRIVATE
    tnt::${PROJECT_NAME}
    Catch2::Catch2WithMain
    fftw::fftw
)

include(Catch)
include(Ctest)

catch_discover_tests(${PROJECT_NAME}_test)