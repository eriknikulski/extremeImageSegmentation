cmake_minimum_required(VERSION 3.21)
project(extremeImageSegmentation C)

set(CMAKE_C_STANDARD 11)


file(GLOB SOURCES
        *.h
        *.c
        )

add_library(libpng STATIC IMPORTED)
set_target_properties(libpng PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/libpng-1.6.37/libpng.a)

add_executable(extremeImageSegmentation ${SOURCE_FILES} ${SOURCES} ${HEADER_FILES})
target_link_libraries (extremeImageSegmentation libpng)
target_link_libraries (extremeImageSegmentation z)
