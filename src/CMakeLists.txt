include(add_scasm_executable)

add_scasm_executable(orbital orbital.c lua_orbital.c
                     nbody.c verlet_ode_solver.c barnes_hut_tree.c body.c)

