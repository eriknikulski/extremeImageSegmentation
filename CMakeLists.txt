cmake_minimum_required(VERSION 3.21)
project(extremeImageSegmentation C)

set(CMAKE_C_STANDARD 11)


file(GLOB SOURCES
        *.h
        *.c
        )

add_executable(extremeImageSegmentation ${SOURCE_FILES} ${SOURCES})
