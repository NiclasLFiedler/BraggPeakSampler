##----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.16...3.27)
project(braggsampler)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()
find_package(DICOM REQUIRED PATHS )

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
#message(${Geant4_DIR})
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(~/geant4-v11.3.0/examples/extended/medical/DICOM/include)
#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#geant4-v11.3.0/examples/extended/medical/DICOM
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB sources_DICOM ~/geant4-v11.3.0/examples/extended/medical/DICOM/src/*.cc)
file(GLOB headers_DICOM ~/geant4-v11.3.0/examples/extended/medical/DICOM/include/*.hh)


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(braggsampler braggsampler.cc ${sources} ${headers} ${sources_DICOM} ${headers_DICOM})
target_link_libraries(braggsampler ${Geant4_LIBRARIES} DICOM)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B2a. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

set(BRAGGSAMPLER_SCRIPTS
  vis.mac
  vis_wrl.mac
  run_p.mac
  ColourMap.dat
  )

foreach(_script ${BRAGGSAMPLER_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS braggsampler DESTINATION bin)
