include(add_polymec_executable)

add_polymec_executable(orbital orbital.c lua_orbital.c 
                       n_squared.c verlet_ode_solver.c body.c)

