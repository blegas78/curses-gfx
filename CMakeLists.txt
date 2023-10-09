#cmake_minimum_required(VERSION 3.0)

cmake_minimum_required(VERSION 3.1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)


function( new_install LIBNAME )
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

if(APPLE)
	set(CMAKE_MACOSX_RPATH 1)
endif()


project(curses-gfx)
set(LIBNAME cursesgfx)

find_package(PNG REQUIRED)
find_package(ASSIMP)



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

set(c_sources src/cursesGfx3d.cpp
                      src/cursesGfx-drawPolygon.cpp
                      src/cursesGfxRenderPipeline.cpp
                      src/cursesGfx-set.cpp
                      src/cursesGfx.cpp
                      src/cursesGfxHandler.cpp
                      src/cursesClock.cpp
                      src/cursesGfxTexture.cpp
                      src/cursesGfxPngLoader.cpp
                      src/cursesGfxTerminal.cpp
                      ) 
set(c_headers  include/curses-gfx-3d.h
		include/curses-gfx.h
		include/curses-gfx-handler.h
		include/curses-clock.h
        include/curses-gfx-texture.h
        include/curses-gfx-png-loader.h
        include/curses-gfx-types.h
        )
        
if(ASSIMP_FOUND)
set(c_sources ${c_sources} src/cursesGfxLoader.cpp )
set(c_headers ${c_headers} include/curses-gfx-loader.h )
set(DEPENDENCIES ${DEPENDENCIES} ${ASSIMP_LIBRARY})
endif(ASSIMP_FOUND)
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

#find_package(PkgConfig REQUIRED)
#pkg_check_modules(ncurses  REQUIRED ncurses )
#if(ncurses_FOUND)
#    message(" --- ncurses_FOUND")
#    message(" --- ${ncurses_LIBRARIES}")
#endif(ncurses_FOUND)


add_executable(asciilines src/testAsciiLines.cpp)
target_link_libraries(asciilines ncurses ${LIBNAME})

add_executable(chaos_clocks src/chaosClocks.cpp)
target_link_libraries(chaos_clocks ncurses ${LIBNAME})

add_executable(stopwatch src/stopwatch.cpp)
target_link_libraries(stopwatch ncurses ${LIBNAME})


add_executable(cube src/cube.cpp)
target_link_libraries(cube ncurses ${LIBNAME})

add_executable(cube2 src/cube2.cpp)
target_link_libraries(cube2 ncurses ${LIBNAME})

add_executable(cube3 src/cube3.cpp)
target_link_libraries(cube3 ncurses ${LIBNAME})

if(PNG_FOUND)
	message("-- Building examples with libpng!")
	add_executable(cubepng src/cubepng.cpp)
	target_link_libraries(cubepng ncurses ${PNG_LIBRARY} ${LIBNAME})
endif(PNG_FOUND)

if(ASSIMP_FOUND)
    message("-- Building examples with libassimp!")
    add_executable(fileLoader src/fileLoader.cpp)
    target_link_libraries(fileLoader ncurses ${LIBNAME} ${ASSIMP_LIBRARY})
endif(ASSIMP_FOUND)

add_executable(aviz src/aviz.cpp)
target_link_libraries(aviz ncurses ${LIBNAME})

add_executable(test3D src/test3D.cpp)
target_link_libraries(test3D ncurses ${LIBNAME})

add_executable(testFill src/testFill.cpp)
target_link_libraries(testFill ncurses ${LIBNAME})

