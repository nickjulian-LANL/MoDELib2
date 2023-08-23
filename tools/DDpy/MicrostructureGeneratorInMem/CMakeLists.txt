include_directories( ${DDPY_INCLUDE_Dir})
include_directories( ${DDPY_INCLUDE_Dir/MicrostructureGeneratorInMem})
include_directories( Eigen3::Eigen)

target_sources(${PROJECT_NAME} PRIVATE
   ${CMAKE_CURRENT_SOURCE_DIR}/MicrostructureSpecification.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/MicrostructureGeneratorInMem.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/MicrostructureGeneratorBaseInMem.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/PeriodicDipoleGeneratorInMem.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/PeriodicLoopGeneratorInMem.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/PrismaticLoopGeneratorInMem.cpp
   )
