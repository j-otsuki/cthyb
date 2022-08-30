# Build hybqmc executable

set(SOURCES
    hyb_qmc_main.cpp
    hyb_qmc.cpp
    ct_qmc_share.cpp
    mt.cpp
    fft_gsl.cpp
    fft_fftw.cpp
    gtau.cpp
    op_determinant.cpp
    operators.cpp
    pade.cpp
)
add_executable(hybqmc ${SOURCES})

# libraries
set(LIBS
    ${MPI_C_LIBRARIES}
    ${GSL_LIBRARIES}
    ${GSL_CBLASLIBRARIES}
    ${FFTW_LIBRARIES}
)
target_link_libraries(hybqmc ${LIBS})

# includes
set(INCLUDES
    ${PROJECT_SOURCE_DIR}/include
    ${GSL_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${FFTW_INCLUDE_DIRS}
)
target_include_directories(hybqmc PRIVATE ${INCLUDES})

install(TARGETS hybqmc)
