#
# Encapsulates building FFTW as an External Project.
include(CheckCCompilerFlag)


set(msg "ATTENTION: You have enabled the use of fftw,")
set(msg "${msg} this library is distributed under a GPL license.")
set(msg "${msg} By enabling this option, the binary of the ITK libraries")
set(msg "${msg} that you are going to build will be covered by a GPL license,")
set(msg "${msg} and so it will be any executable that you link against these libraries.")
message("${msg}")

#--check_c_compiler_flag(-fopenmp C_HAS_fopenmp)
#--if(${C_HAS_fopenmp} AND FALSE)
#--    set(FFTW_THREADS_CONFIGURATION --enable-openmp)
#--    set(OPENMP_FLAG "-fopenmp")
#--  else()
    set(FFTW_THREADS_CONFIGURATION --enable-threads)
    set(OPENMP_FLAG "")
#--endif()

#--Some influential environment variables:
#--  CC          C compiler command
#--  CFLAGS      C compiler flags
#--  LDFLAGS     linker flags, e.g. -L<lib dir> if you have libraries in a
#--              nonstandard directory <lib dir>
#--  LIBS        libraries to pass to the linker, e.g. -l<library>
#--  CPPFLAGS    C/C++/Objective C preprocessor flags, e.g. -I<include dir> if
#--              you have headers in a nonstandard directory <include dir>
#-- set(ENV{CC}       "${CMAKE_C_COMPILER}")
#-- set(ENV{CFLAGS}   "${CMAKE_C_FLAGS} ${OPENMP_FLAG}")
#-- set(ENV{LDFLAGS}  "${CMAKE_C_FLAGS} ${OPENMP_FLAG}")
#-- set(ENV{LIBS}     "${CMAKE_EXE_LINKER_FLAGS} ${OPENMP_FLAG}")
#-- set(ENV{CPPFLAGS} "${CMAKE_C_FLAGS} ${OPENMP_FLAG}")

## Perhaps in the future a set of TryCompiles could be used here.
set(FFTW_OPTIMIZATION_CONFIGURATION "" CACHE INTERNAL "architecture flags: --enable-sse --enable-sse2 --enable-altivec --enable-mips-ps --enable-cell")
if(USE_SYSTEM_FFTW)
  find_package( FFTW )
  link_directories(${FFTW_LIBDIR})
else(USE_SYSTEM_FFTW)
  if(WIN32 AND NOT MINGW)
    message("Can't build fftw as external project on Windows")
    message(ERROR "install fftw and use USE_SYSTEM_FFTW")
  else(WIN32 AND NOT MINGW)
    #
    # fftw limitation -- can't be built in
    # a directory with whitespace in its name.
    if(${CMAKE_CURRENT_BINARY_DIR} MATCHES ".*[ \t].*")
      message(FATAL_ERROR
        "Can't build fftw in a directory with whitespace in its name")
    endif(${CMAKE_CURRENT_BINARY_DIR} MATCHES ".*[ \t].*")
    #
    # build fftw as an external project
    if(BUILD_SHARED_LIBS)
      set(FFTW_SHARED_FLAG --enable-shared)
    endif(BUILD_SHARED_LIBS)
    if(USE_FFTWF)
      ExternalProject_add(fftwf
        PREFIX fftwf
        URL "http://www.fftw.org/fftw-3.2.2.tar.gz"
        URL_MD5 b616e5c91218cc778b5aa735fefb61ae
        CONFIGURE_COMMAND ${ITK_BINARY_DIR}/fftwf/src/fftwf/configure
        ${FFTW_SHARED_FLAG}
        ${FFTW_OPTIMIZATION_CONFIGURATION}
        ${FFTW_THREADS_CONFIGURATION}
        --disable-fortran
        --enable-float
        --prefix=${ITK_BINARY_DIR}/fftw
        )
    endif(USE_FFTWF)

    if(USE_FFTWD)
      ExternalProject_add(fftwd
        PREFIX fftwd
        URL "http://www.fftw.org/fftw-3.2.2.tar.gz"
        URL_MD5 b616e5c91218cc778b5aa735fefb61ae
        CONFIGURE_COMMAND ${ITK_BINARY_DIR}/fftwd/src/fftwd/configure
        ${FFTW_SHARED_FLAG}
        ${FFTW_OPTIMIZATION_CONFIGURATION}
        ${FFTW_THREADS_CONFIGURATION}
        --disable-fortran
        --disable-float
        --prefix=${ITK_BINARY_DIR}/fftw
        )
    endif(USE_FFTWD)
    link_directories(${ITK_BINARY_DIR}/fftw/lib)
    include_directories(${ITK_BINARY_DIR}/fftw/include)
    # backwards compatibility
    set(FFTW_INCLUDE_PATH ${ITK_BINARY_DIR}/fftw/include)
    #
    # copy libraries into install tree
    install(CODE
      "file(GLOB FFTW_LIBS ${ITK_BINARY_DIR}/fftw/lib/*fftw3*)
file(INSTALL DESTINATION \"\${CMAKE_INSTALL_PREFIX}/lib/InsightToolkit-${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}\"
TYPE FILE FILES \${FFTW_LIBS})")
    #
    # copy headers into install tree
    install(CODE
      "file(GLOB FFTW_INC ${ITK_BINARY_DIR}/fftw/include/*fftw3*)
file(INSTALL DESTINATION \"\${CMAKE_INSTALL_PREFIX}/include/InsightToolkit-${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}/Algorithms\"
TYPE FILE FILES \${FFTW_INC})")

  endif(WIN32 AND NOT MINGW)
endif(USE_SYSTEM_FFTW)
