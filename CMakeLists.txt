cmake_minimum_required(VERSION 3.5)
project(px_codec)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    src/bitio/bit.h
    src/bitio/bitarray.h
    src/bitio/bitdebug.h
    src/bitio/bitfile.h
    src/bitio/bitview.h
    src/compression/arithmetic_coder.h
    src/compression/bitmanip.h
    src/compression/compression.h
    src/compression/elias_gamma.h
    src/compression/huffman.h
    src/compression/math.h
    src/compression/rice.h
    src/compression/runlength.h
    src/lodepng.cpp
    src/lodepng.h
    src/main.cpp
    src/util.h)

add_executable(px_codec ${SOURCE_FILES})

IF (WIN32)
    TARGET_LINK_LIBRARIES(px_codec)
ELSEIF (UNIX)
    TARGET_LINK_LIBRARIES(px_codec
            gsl
            gslcblas)
ENDIF ()