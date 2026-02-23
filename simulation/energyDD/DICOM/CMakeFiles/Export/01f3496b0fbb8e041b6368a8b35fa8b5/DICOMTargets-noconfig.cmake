#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "DICOM-library" for configuration ""
set_property(TARGET DICOM-library APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(DICOM-library PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libDICOM.so"
  IMPORTED_SONAME_NOCONFIG "libDICOM.so"
  )

list(APPEND _cmake_import_check_targets DICOM-library )
list(APPEND _cmake_import_check_files_for_DICOM-library "${_IMPORT_PREFIX}/lib/libDICOM.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
