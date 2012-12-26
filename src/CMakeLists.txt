set(TOV_solver_SOURCES
    eos_pres.c
    eos_rho.c
    func.c
    main.c
    make_grid.c
    root_func.c
    )

add_executable(TOV_solver ${TOV_solver_SOURCES})

target_link_libraries(TOV_solver
    m
    gsl
    gslcblas
    )
