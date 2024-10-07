
cmake_minimum_required (VERSION 3.12)

# double conversion
add_subdirectory (double-conversion)
set (DOUBLE_CONVERSION_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/double-conversion PARENT_SCOPE)
set_property (TARGET double-conversion PROPERTY FOLDER "thirdparty")

# glm
set (GLM_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/glm PARENT_SCOPE)