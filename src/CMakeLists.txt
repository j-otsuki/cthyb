
set(SOURCES
    hyb_qmc_main.cpp
    hyb_qmc.cpp
    ct_qmc_share.cpp
    mt.cpp
    fft.cpp
    gtau.cpp
    op_determinant.cpp
    operators.cpp
    pade.cpp
)
add_executable(hybqmc ${SOURCES})

set(LIBS
    ${MPI_C_LIBRARIES}
    ${GSL_LIBRARIES}
    ${GSL_CBLASLIBRARIES}
)
target_link_libraries(hybqmc ${LIBS})

target_include_directories(hybqmc PRIVATE ${PROJECT_SOURCE_DIR}/include)

install(TARGETS hybqmc)
