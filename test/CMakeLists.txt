Include(FetchContent)

FetchContent_Declare(Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.4.0
)

FetchContent_MakeAvailable(Catch2)

include(CTest)
include(Catch)

add_executable(${PROJECT_NAME}_test
    main.cpp
    fourier_transform.cpp
    
    
)

target_compile_options(${PROJECT_NAME}_test PRIVATE -D_USE_MATH_DEFINES)

target_link_libraries(${PROJECT_NAME}_test PRIVATE
    tnt::${PROJECT_NAME}
    Catch2::Catch2WithMain
)

catch_discover_tests(${PROJECT_NAME}_test)