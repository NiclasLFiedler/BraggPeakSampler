set(DICOM_VERSION 11.3.0)

macro(set_and_check _var _file)
    set(${_var} "${_file}")
    if(NOT EXISTS "${_file}")
        message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
    endif()
endmacro()

macro(check_required_components _NAME)
    foreach(comp ${${_NAME}_FIND_COMPONENTS})
        if(NOT ${_NAME}_${comp}_FOUND)
            if(${_NAME}_FIND_REQUIRED_${comp})
                set(${_NAME}_FOUND FALSE)
            endif()
        endif()
    endforeach()
endmacro()


set_and_check(DICOM_INCLUDE_DIR "/home/niclas/geant4-v11.3.0/examples/extended/medical/DICOM/include")

check_required_components(DICOM)

include(${CMAKE_CURRENT_LIST_DIR}/DICOMBuild.cmake)

set(DICOM_INCLUDE_DIRS ${DICOM_INCLUDE_DIR})
set(DICOM_LIBRARIES DICOM-library)

set(DICOM_USE_DCMTK OFF)
set(DICOM_USE_HEAD OFF)

if(DICOM_USE_DCMTK)
    set(dicomReader_DIR ${CMAKE_CURRENT_LIST_DIR}
        CACHE PATH "Path to dicomReader configuration")
    find_package(dicomReader 11.3.0 EXACT REQUIRED)
    list(APPEND DICOM_INCLUDE_DIRS ${dicomReader_INCLUDE_DIRS})
    list(APPEND DICOM_LIBRARIES ${dicomReader_LIBRARIES})
    add_definitions(-DG4_DCMTK)
endif(DICOM_USE_DCMTK)

if(DICOM_USE_HEAD)
    add_definitions(-DDICOM_USE_HEAD)
endif(DICOM_USE_HEAD)

set(DICOM_FOUND ON)
