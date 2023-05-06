CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
PROJECT(fit3d) #project name

######################################################
##################ADJUST##############################
######################################################
#what shall be build?
SET(MEXF TRUE) #generate mex file for matlab usage
SET(CUDA TRUE) #use cuda for speedup if applicable

##################CPP & CUDA FILES####################
######################################################

#cuda src files (e.g. kernels)
file(GLOB_RECURSE CU_SRC_FILES "GPUmleFit_LM/src/*" )

#cuda lib files (e.g. kernels)
SET(CU_LIB_FILES "GPUmleFit_LM/lib/" )

#cpp lib files 
SET(CPP_LIB_FILES "CPUmleFit_LM/lib/")

#cpp src files 
file(GLOB_RECURSE CPP_SRC_FILES "CPUmleFit_LM/src/*")

######################################################
##################SET EVERYTHING UP###################
######################################################

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)    #output for mex-file i.e., MEXF=TRUE

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake) #add modules e.g., FindMatlab.cmake

MESSAGE(STATUS "MEXF=" ${MEXF})
MESSAGE(STATUS "CUDA=" ${CUDA})


##################FIND MATLAB#########################
######################################################
if(${MEXF})
  #not sure if the following two definitions are necessary
  ADD_DEFINITIONS(/DMATLAB_MEX_FILE) #define matlab macros
  ADD_DEFINITIONS(/DMX_COMPAT_32)

  FIND_PACKAGE(Matlab REQUIRED)
  #add definition USE_MEX for C++ code, to be able to use preprocessor if matlab is not used
  ADD_DEFINITIONS(-DUSE_MEX)

  # set up matlab libraries
  INCLUDE_DIRECTORIES(${Matlab_INCLUDE_DIRS})

  #use correct suffix depending on OS
  if(WIN32) # 32-bit or 64-bit mex
      if (CMAKE_CL_64)
        SET( MEX_SUFFIX .mexw64 )
      else(CMAKE_CL_64)
        SET( MEX_SUFFIX .mexw32 )
      endif(CMAKE_CL_64)
  else(WIN32)
      if (CMAKE_SIZEOF_VOID_P MATCHES "8")
        SET( MEX_SUFFIX .mexa64 )
      else(CMAKE_SIZEOF_VOID_P MATCHES "8")
        SET( MEX_SUFFIX .mexglx )
      endif (CMAKE_SIZEOF_VOID_P MATCHES "8")
  endif(WIN32)
endif(${MEXF})

##################FIND CUDA###########################
######################################################
if(${CUDA})
  FIND_PACKAGE(CUDA REQUIRED)
  if(CUDA_FOUND)  
    message(STATUS "CUDA Found, gpu version will be compiled.")
    #add definition USE_MEX for C++ code, to be able to use preprocessor if CUDA is not used
    ADD_DEFINITIONS(-DUSE_CUDA)
    include_directories(${CUDA_INCLUDE_DIRS} ${CU_LIB_FILES})
    
    #file(s) for mex ( e.g. mex entrance function) version only
    SET(MEX_MAIN_FILES GPUmleFit_LM/mexFunction.cpp)


    # set SRC_FILES using the common cpp files and cu source files
    SET(SRC_FILES ${CU_SRC_FILES})
  else(CUDA_FOUND)
    #cpu version
    message("CUDA not Found, cpu version will be compiled.")
    SET(CUDA FALSE)

  endif(CUDA_FOUND)
endif(${CUDA})
if(NOT ${CUDA})
  #cpu version
  message("-- cpu version will be compiled.")
  #file(s) for mex ( e.g. mex entrance function) version only
  SET(MEX_MAIN_FILES CPUmleFit_LM/mexFunction.cpp)
  include_directories(${CPP_LIB_FILES})
  SET(SRC_FILES ${CPP_SRC_FILES})
endif(NOT ${CUDA})

######################################################
##################BUILD FILES#########################
######################################################

######################MEX#############################
######################################################
if(${MEXF})
 
  if(${CUDA})
    SET(MEX_LIB_NAME "${PROJECT_NAME}GPU")
    CUDA_COMPILE(${MEX_LIB_NAME} ${MEX_MAIN_FILES} ${SRC_FILES} SHARED)   
  else(${CUDA})
    SET(MEX_LIB_NAME "${PROJECT_NAME}CPU")
  endif(${CUDA})

  
  ADD_LIBRARY(${MEX_LIB_NAME} SHARED ${${MEX_LIB_NAME}} ${MEX_MAIN_FILES} ${SRC_FILES} ${CMAKE_SOURCE_DIR}/matlab.def)
  message(${MEX_LIB_NAME})
  if(${CUDA})
    TARGET_LINK_LIBRARIES(${MEX_LIB_NAME}
    			   ${Matlab_LIBRARIES}
                          ${CUDA_LIBRARIES} 
                          )
    SET_TARGET_PROPERTIES(${MEX_LIB_NAME} PROPERTIES LINKER_LANGUAGE CXX)
  else(${CUDA})
  TARGET_LINK_LIBRARIES(${MEX_LIB_NAME}
                        ${Matlab_LIBRARIES}
                        )
  endif(${CUDA})
    SET_TARGET_PROPERTIES(${MEX_LIB_NAME} PROPERTIES PREFIX "" SUFFIX ${MEX_SUFFIX})
endif(${MEXF})

