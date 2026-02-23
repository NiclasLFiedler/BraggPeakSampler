# Install script for directory: /home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/DICOM" TYPE FILE FILES
    "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/include/DicomActionInitialization.hh"
    "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/include/DicomDetectorConstruction.hh"
    "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/include/DicomEventAction.hh"
    "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/include/DicomHandler.hh"
    "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/include/DicomIntersectVolume.hh"
    "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/include/DicomNestedParamDetectorConstruction.hh"
    "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/include/DicomNestedPhantomParameterisation.hh"
    "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/include/DicomPartialDetectorConstruction.hh"
    "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/include/DicomPhantomParameterisationColour.hh"
    "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/include/DicomPhantomZSliceHeader.hh"
    "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/include/DicomPhantomZSliceMerged.hh"
    "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/include/DicomPrimaryGeneratorAction.hh"
    "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/include/DicomRegularDetectorConstruction.hh"
    "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/include/DicomRun.hh"
    "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/include/DicomRunAction.hh"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DICOM" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DICOM")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DICOM"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/niclas/braggpeak_sampler/simulation/detector/DICOM/DICOM")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DICOM" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DICOM")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DICOM"
         OLD_RPATH "/home/niclas/braggpeak_sampler/simulation/detector/DICOM:/usr/local/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DICOM")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/niclas/braggpeak_sampler/simulation/detector/DICOM/CMakeFiles/DICOM.dir/install-cxx-module-bmi-noconfig.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDICOM.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDICOM.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDICOM.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/niclas/braggpeak_sampler/simulation/detector/DICOM/libDICOM.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDICOM.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDICOM.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDICOM.so"
         OLD_RPATH "/usr/local/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDICOM.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/niclas/braggpeak_sampler/simulation/detector/DICOM/CMakeFiles/DICOM-library.dir/install-cxx-module-bmi-noconfig.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/build/DICOMTargets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/build/DICOMTargets.cmake"
         "/home/niclas/braggpeak_sampler/simulation/detector/DICOM/CMakeFiles/Export/01f3496b0fbb8e041b6368a8b35fa8b5/DICOMTargets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/build/DICOMTargets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/build/DICOMTargets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/build/DICOMTargets.cmake")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/build" TYPE FILE FILES "/home/niclas/braggpeak_sampler/simulation/detector/DICOM/CMakeFiles/Export/01f3496b0fbb8e041b6368a8b35fa8b5/DICOMTargets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^()$")
    list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
     "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/build/DICOMTargets-noconfig.cmake")
    if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
      message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
    endif()
    if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
      message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
    endif()
    file(INSTALL DESTINATION "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/build" TYPE FILE FILES "/home/niclas/braggpeak_sampler/simulation/detector/DICOM/CMakeFiles/Export/01f3496b0fbb8e041b6368a8b35fa8b5/DICOMTargets-noconfig.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/build/DICOMConfig.cmake;/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/build/DICOMConfigVersion.cmake")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/build" TYPE FILE FILES
    "/home/niclas/braggpeak_sampler/simulation/detector/DICOM/InstallTreeFiles/DICOMConfig.cmake"
    "/home/niclas/braggpeak_sampler/simulation/detector/DICOM/InstallTreeFiles/DICOMConfigVersion.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/niclas/braggpeak_sampler/simulation/detector/DICOM/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
