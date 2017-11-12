// Copyright (c) 2016, Jeffrey N. Johnson
// All rights reserved.

#include "constants.h"
#include "lua_orbital.h"

#include "lualib.h"
#include "lauxlib.h"

static int b_new(lua_State* L)
{
  if ((lua_gettop(L) != 1) || !lua_istable(L, 1))
    luaL_error(L, "Argument must be a table containing name, x, v.");

  lua_getfield(L, 1, "name");
  if (!lua_isstring(L, -1))
    luaL_error(L, "name must be a string.");
  const char* name = lua_tostring(L, -1);

  lua_getfield(L, 1, "m");
  if (!lua_is_real(L, -1))
    luaL_error(L, "m must be a positive number.");
  real_t m = lua_to_real(L, -1);
  if (m <= 0.0)
    luaL_error(L, "m must be positive.");

  lua_getfield(L, 1, "x");
  if (!lua_is_point(L, -1))
    luaL_error(L, "x must be a point.");
  point_t* x = lua_to_point(L, -1);

  lua_getfield(L, 1, "v");
  if (!lua_is_vector(L, -1))
    luaL_error(L, "v must be a vector.");
  vector_t* v = lua_to_vector(L, -1);

  body_t* b = body_new(name, m, x, v);
  lua_push_body(L, b);
  return 1;
}

static lua_module_function body_functions[] = {
  {"new", b_new, "Creates a new massive body with name, m, x, v."},
  {NULL, NULL, NULL}
};

static int b_name(lua_State* L)
{
  body_t* b = lua_to_body(L, 1);
  lua_pushstring(L, b->name);
  return 1;
}

static int b_mass(lua_State* L)
{
  body_t* b = lua_to_body(L, 1);
  lua_push_real(L, b->m);
  return 1;
}

static int b_x(lua_State* L)
{
  body_t* b = lua_to_body(L, 1);
  lua_push_point(L, point_new(b->x.x, b->x.y, b->x.z));
  return 1;
}

static int b_v(lua_State* L)
{
  body_t* b = lua_to_body(L, 1);
  lua_push_vector(L, vector_new(b->v.x, b->v.y, b->v.z));
  return 1;
}

static lua_record_field body_fields[] = {
  {"name", b_name, NULL},
  {"mass", b_mass, NULL},
  {"x", b_x, NULL},
  {"v", b_v, NULL},
  {NULL, NULL, NULL}
};

static int body_tostring(lua_State* L)
{
  body_t* b = lua_to_body(L, 1);
  lua_pushfstring(L, "body (name = %s, m = %f, x = (%f, %f, %f), v = (%f, %f, %f))", 
                  b->name, b->m, b->x.x, b->x.y, b->x.z, b->v.x, b->v.y, b->v.z);
  return 1;
}

static lua_record_metamethod body_mm[] = {
  {"__tostring", body_tostring},
  {NULL, NULL}
};

static void get_nbody_args(lua_State* L, real_t* G, body_array_t** bodies)
{
  // Where all our bodies at?
  lua_getfield(L, 1, "bodies");
  if (!lua_istable(L, -1))
    luaL_error(L, "bodies must be a list of body objects.");
  *bodies = body_array_new();
  {
    int i = 1;
    while (true)
    {
      lua_rawgeti(L, -1, i);
      if (lua_isnil(L, -1)) break;
      if (!lua_is_body(L, -1))
        luaL_error(L, "Item %d in bodies is not a body.", i);
      body_t* b = lua_to_body(L, -1);
      body_array_append_with_dtor(*bodies, b, body_free);
      lua_pop(L, 1);
      ++i;
    }
  }
  lua_pop(L, 1);

  // Are we given a gravitational constant G?
  *G = GRAVITATIONAL_CONSTANT;
  lua_getfield(L, 1, "G");
  if (lua_isnumber(L, -1))
  {
    *G = lua_to_real(L, -1);
    if (*G <= 0.0)
      luaL_error(L, "G must be positive.");
  }
}

// This function constructs our brute-force n-body model. It takes a table 
// of arguments.
static int nb_brute_force(lua_State* L)
{
  if ((lua_gettop(L) != 1) || !lua_istable(L, 1))
    luaL_error(L, "Argument must be a table of parameters and bodies.");

  // Fetch our parameters.
  real_t G;
  body_array_t* bodies;
  get_nbody_args(L, &G, &bodies);

  // Create the thingy and push it to the stack.
  model_t* nb = brute_force_nbody_new(G, bodies);

  lua_push_model(L, nb);
  return 1;
}

// This function constructs our Barnes-Hut hierarchical n-body model. It 
// takes a table of arguments.
static int nb_barnes_hut(lua_State* L)
{
  if ((lua_gettop(L) != 1) || !lua_istable(L, 1))
    luaL_error(L, "Argument must be a table of parameters and bodies.");

  // Fetch our parameters.
  real_t G;
  body_array_t* bodies;
  get_nbody_args(L, &G, &bodies);

  // Are we given a far-field parameter theta?
  real_t theta = 0.5;
  lua_getfield(L, 1, "theta");
  if (lua_isnumber(L, -1))
  {
    theta = lua_to_real(L, -1);
    if (theta < 0.0)
      luaL_error(L, "theta must be non-negative.");
  }

  // Create the thingy and push it to the stack.
  model_t* nb = barnes_hut_nbody_new(G, theta, bodies);

  lua_push_model(L, nb);
  return 1;
}

static int nb_x_probe(lua_State* L)
{
  if (!lua_isstring(L, 1))
    luaL_error(L, "Argument must be the name of a body.");
  const char* name = lua_tostring(L, 1);
  lua_push_probe(L, nbody_x_probe_new(name));
  return 1;
}

static int nb_v_probe(lua_State* L)
{
  if (!lua_isstring(L, 1))
    luaL_error(L, "Argument must be the name of a body.");
  const char* name = lua_tostring(L, 1);
  lua_push_probe(L, nbody_v_probe_new(name));
  return 1;
}

static int nb_E_probe(lua_State* L)
{
  lua_push_probe(L, nbody_E_probe_new());
  return 1;
}

static void lua_register_nbody(lua_State* L)
{
  // Create a new table and fill it with our gravity models.
  lua_newtable(L);

  lua_pushcfunction(L, nb_brute_force);
  lua_setfield(L, -2, "brute_force");

  lua_pushcfunction(L, nb_barnes_hut);
  lua_setfield(L, -2, "barnes_hut");

  lua_pushcfunction(L, nb_x_probe);
  lua_setfield(L, -2, "x_probe");

  lua_pushcfunction(L, nb_v_probe);
  lua_setfield(L, -2, "v_probe");

  lua_pushcfunction(L, nb_E_probe);
  lua_setfield(L, -2, "E_probe");

  lua_setglobal(L, "nbody");
}

static void lua_register_constants(lua_State* L)
{
  // Create a new table and fill it with our named physical constants.
  lua_newtable(L);

  lua_pushnumber(L, GRAVITATIONAL_CONSTANT);
  lua_setfield(L, -2, "G");
  lua_pushnumber(L, SOLAR_MASS);
  lua_setfield(L, -2, "solar_mass");
  lua_pushnumber(L, EARTH_MASS);
  lua_setfield(L, -2, "earth_mass");
  lua_pushnumber(L, SPEED_OF_LIGHT);
  lua_setfield(L, -2, "c");
  lua_pushnumber(L, ASTRONOMICAL_UNIT);
  lua_setfield(L, -2, "AU");

  lua_setglobal(L, "constants");
}


int lua_register_orbital_modules(lua_State* L)
{
  lua_register_core_modules(L);
  lua_register_model_modules(L);

  lua_register_record_type(L, "body", "A body in 3D space.", body_functions, body_fields, body_mm);

  lua_register_nbody(L);
  lua_register_constants(L);
  return 0;
}

void lua_push_body(lua_State* L, body_t* body)
{
  lua_push_record(L, "body", body, NULL);
}

bool lua_is_body(lua_State* L, int index)
{
  return lua_is_record(L, index, "body");
}

body_t* lua_to_body(lua_State* L, int index)
{
  return (body_t*)lua_to_record(L, index, "body");
}

