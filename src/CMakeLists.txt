tnt_project__add_library(TARGET ${PROJECT_NAME}
	impl/utility.cpp
)

target_compile_features(${PROJECT_NAME}	PUBLIC cxx_std_20)

target_compile_options(${PROJECT_NAME} PUBLIC -D_USE_MATH_DEFINES)