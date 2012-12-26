set(TOV_solver_SOURCES
    main.c
    )

add_executable(TOV_solver ${TOV_solver_SOURCES})

target_link_libraries(TOV_solver
    m
    gsl
    gslcblas
    )
