cmake_minimum_required(VERSION 3.0)
project(asciilines)

set(CMAKE_BUILD_TYPE Debug)

include_directories(include)

include(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX("linux/fb.h" FB_EXISTS)
if(FB_EXISTS)
#	message("adding definition FB_SUPPORT")
#	add_definitions(-DFB_SUPPORT)
endif(FB_EXISTS)

add_library(cursesgfx src/cursesGfx3d.cpp
                      src/cursesGfx-drawPolygon.cpp
                      src/cursesGfxRenderPipeline.cpp
                      src/cursesGfx-set.cpp
                      include/curses-gfx-3d.h
                      src/cursesGfx.cpp
                      include/curses-gfx.h
                      src/cursesGfxHandler.cpp
                      include/curses-gfx-handler.h
                      src/cursesClock.cpp
                      include/curses-clock.h )
target_link_libraries(cursesgfx ncurses )


add_executable(asciilines src/testAsciiLines.cpp)
target_link_libraries(asciilines ncurses cursesgfx)

add_executable(chaos_clocks src/chaosClocks.cpp)
target_link_libraries(chaos_clocks ncurses cursesgfx)

add_executable(stopwatch src/stopwatch.cpp)
target_link_libraries(stopwatch ncurses cursesgfx)


add_executable(cube src/cube.cpp)
target_link_libraries(cube ncurses cursesgfx)

add_executable(cube2 src/cube2.cpp)
target_link_libraries(cube2 ncurses cursesgfx)

add_executable(test3D src/test3D.cpp)
target_link_libraries(test3D ncurses cursesgfx)

