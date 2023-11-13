#cmake_minimum_required(VERSION 3.0)

cmake_minimum_required(VERSION 3.1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)


function( new_install LIBNAME )
    # Resources
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/resources")
        file(GLOB ${LIBNAME}_RESOURCES resources/*)
        set(${LIBNAME}_resources ${${LIBNAME}_RESOURCES})
        set(resources_subdirectory etc/${LIBNAME} PARENT_SCOPE)
        set(resources_subdirectory etc/${LIBNAME})
        file(COPY
                    ${CMAKE_CURRENT_SOURCE_DIR}/resources
                    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/
                    )
        install( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/resources #${${LIBNAME}_resources}
                    DESTINATION ${resources_subdirectory}
                    COMPONENT resources
                        )
        add_definitions(-DRESOURCES_DIRECTORY="${CMAKE_INSTALL_PREFIX}/${resources_subdirectory}")
        set(CURSES_ROUSOURCE_DIR "${CMAKE_INSTALL_PREFIX}/${resources_subdirectory}/resources/")
        configure_file(include/curses-gfx-resources.h.in "${CMAKE_CURRENT_BINARY_DIR}/curses-gfx-resources.h" @ONLY)
    endif()
    
	# Install
	if(BUILD_SHARED_LIBS)
		install( TARGETS ${LIBNAME}
			LIBRARY DESTINATION lib #lib${LIB_SUFFIX}
			COMPONENT sharedlibs
			)
	endif(BUILD_SHARED_LIBS)

	if(BUILD_STATIC_LIBS)
		install( TARGETS ${LIBNAME}-static
			ARCHIVE DESTINATION lib #lib${LIB_SUFFIX}
			COMPONENT staticlibs
			)
	endif(BUILD_STATIC_LIBS)

	install( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ #${c_headers}
		DESTINATION include
		COMPONENT headers
		)
		
	install( FILES ${GENERATED_HEADERS}
		DESTINATION include
		COMPONENT headers
		)

	message("-- Install Directory ${LIBNAME}: lib${LIB_SUFFIX}")
	message("-- Install Directory ${LIBNAME}-static: lib${LIB_SUFFIX}")
	message("-- Install Directory ${c_headers}: include/${LIB_SUFFIX}")
endfunction()


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_COLOR_MAKEFILE ON)

set(CMAKE_SKIP_BUILD_RPATH false)
set(CMAKE_BUILD_WITH_INSTALL_RPATH false)
#SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib64")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)

if(APPLE)
	set(CMAKE_MACOSX_RPATH 1)
endif()


project(curses-gfx)
set(LIBNAME cursesgfx)

find_package(PNG REQUIRED)
find_package(ASSIMP)
#find_package(BULLET)
find_package(PkgConfig REQUIRED)
pkg_check_modules(BULLET bullet)
pkg_check_modules(SDL2 SDL2)


set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Create suffix to eventually install in lib64
IF(CMAKE_SIZEOF_VOID_P EQUAL 4)
    SET(LIB_SUFFIX "")
    SET(PACK_ARCH "")
  ELSE(CMAKE_SIZEOF_VOID_P EQUAL 8)
    SET(LIB_SUFFIX 64)
    SET(PACK_ARCH .x86_64)
endif(CMAKE_SIZEOF_VOID_P EQUAL 4)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")#${LIB_SUFFIX}")

option( BUILD_SHARED_LIBS "Build the shared library" ON)
option( BUILD_STATIC_LIBS "Build the static library" ON)

find_package(Threads REQUIRED)

include_directories(include)
if(PNG_FOUND)
#    message(" ----- ${PNG_INCLUDE_DIR}")
	include_directories(${PNG_INCLUDE_DIR})
endif(PNG_FOUND)


include(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX("linux/fb.h" FB_EXISTS)
if(FB_EXISTS)
#	message("adding definition FB_SUPPORT")
#	add_definitions(-DFB_SUPPORT)
endif(FB_EXISTS)

if(EXISTS /home/pi/rpi-rgb-led-matrix)
	message("adding definition RGB_MATRIX_SUPPORT")
	add_definitions(-DRGB_MATRIX_SUPPORT)
	
endif()

set(DEPENDENCIES ${DEPENDENCIES} ncurses ${CMAKE_THREAD_LIBS_INIT} ${PNG_LIBRARY})


set(GENERATED_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/curses-gfx-resources.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(c_sources src/cursesGfx3d.cpp
                      src/cursesGfx-drawPolygon.cpp
                      src/cursesGfxRenderPipeline.cpp
                      src/cursesGfx-set.cpp
                      src/cursesGfx.cpp
                      src/cursesGfxTypes.cpp
                      src/cursesGfxHandler.cpp
                      src/cursesClock.cpp
                      src/cursesGfxTexture.cpp
                      src/cursesGfxPngLoader.cpp
                      src/cursesGfxTerminal.cpp
                      src/cursesGfxThreadpool.cpp
                      ) 
set(c_headers  include/curses-gfx-3d.h
		include/curses-gfx.h
		include/curses-gfx-handler.h
		include/curses-clock.h
        include/curses-gfx-texture.h
        include/curses-gfx-png-loader.h
        include/curses-gfx-types.h
        ${GENERATED_HEADERS}
        )
        
if(ASSIMP_FOUND)
    set(c_sources ${c_sources} src/cursesGfxLoader.cpp )
    set(c_headers ${c_headers} include/curses-gfx-loader.h )
    set(DEPENDENCIES ${DEPENDENCIES} ${ASSIMP_LIBRARY})
endif(ASSIMP_FOUND)

if(BULLET_FOUND)
    message(" - BULLET_INCLUDE_DIRS: ${BULLET_INCLUDE_DIRS}")
    include_directories(${BULLET_INCLUDEDIR})
    set(c_sources ${c_sources} src/cursesGfxPhysics.cpp )
    set(c_headers ${c_headers} include/curses-gfx-physics.h )
    link_directories(${BULLET_LIBDIR})
    set(DEPENDENCIES ${DEPENDENCIES} ${BULLET_LIBRARIES})
endif(BULLET_FOUND)

if(SDL2_FOUND)
    link_directories(${SDL2_LIBDIR})
endif(SDL2_FOUND)
set(c_sources ${c_sources} ${c_headers} )

# Libraries:
#add_library(${LIBNAME} ${c_sources} )

add_library(${LIBNAME}  SHARED   ${c_sources} )
add_library(${LIBNAME}-static STATIC ${c_sources})
set_target_properties(${LIBNAME}-static PROPERTIES OUTPUT_NAME "${LIBNAME}")
set_target_properties(${LIBNAME} PROPERTIES CLEAN_DIRECT_OUTPUT 1)
set_target_properties(${LIBNAME}-static PROPERTIES CLEAN_DIRECT_OUTPUT 1)

target_link_libraries(${LIBNAME} ${DEPENDENCIES})


new_install( ${LIBNAME} )



add_executable(asciilines examples/testAsciiLines.cpp)
target_link_libraries(asciilines ncurses ${LIBNAME})

add_executable(chaos_clocks examples/chaosClocks.cpp)
target_link_libraries(chaos_clocks ncurses ${LIBNAME})

add_executable(stopwatch examples/stopwatch.cpp)
target_link_libraries(stopwatch ncurses ${LIBNAME})


add_executable(cube examples/cube.cpp)
target_link_libraries(cube ncurses ${LIBNAME})

add_executable(cube2 examples/cube2.cpp)
target_link_libraries(cube2 ncurses ${LIBNAME})

add_executable(cube3 examples/cube3.cpp)
target_link_libraries(cube3 ncurses ${LIBNAME})

if(PNG_FOUND)
	message("-- Building examples with libpng!")
	add_executable(cubepng examples/cubepng.cpp)
	target_link_libraries(cubepng ncurses ${PNG_LIBRARY} ${LIBNAME})
endif(PNG_FOUND)

if(ASSIMP_FOUND)
    message("-- Building examples with libassimp!")
    add_executable(fileLoader examples/fileLoader.cpp)
    target_link_libraries(fileLoader ncurses ${LIBNAME} ${ASSIMP_LIBRARY})
endif(ASSIMP_FOUND)

if(BULLET_FOUND)
    message("-- Building examples with bullet physics!")
    add_executable(bullet examples/bullet.cpp)
    target_link_libraries(bullet ncurses ${LIBNAME} ${BULLET_LIBRARY})
endif(BULLET_FOUND)

if(ASSIMP_FOUND AND BULLET_FOUND AND SDL2_FOUND)
    add_executable(game examples/game.cpp)
#    message( " -------- SDL2_LIBRARY = ${SDL2_LIBRARIES}")
    target_link_libraries(game ncurses ${LIBNAME} ${ASSIMP_LIBRARY} ${BULLET_LIBRARY} ${SDL2_LIBRARIES})
endif(ASSIMP_FOUND AND BULLET_FOUND AND SDL2_FOUND)

add_executable(aviz examples/aviz.cpp)
target_link_libraries(aviz ncurses ${LIBNAME})

add_executable(testColor examples/testColor.cpp)
target_link_libraries(testColor ncurses ${LIBNAME})

add_executable(testFill examples/testFill.cpp)
target_link_libraries(testFill ncurses ${LIBNAME})


add_executable(threadPoolCheck examples/threadPoolCheck.cpp)
target_link_libraries(threadPoolCheck ncurses ${LIBNAME})

if(ASSIMP_FOUND)
    add_executable(deferred examples/deferred.cpp)
    target_link_libraries(deferred ncurses ${LIBNAME})
endif(ASSIMP_FOUND)
