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

static int b_sc(lua_State* L)
{
  b_new(L);
  body_t* b = lua_to_body(L, -1);
  b->schwartzchild = true;
  return 1;
}

static lua_module_function body_functions[] = {
  {"new", b_new, "Creates a new massive body with name, m, x, v."},
  {"schwartzchild", b_sc, "Creates a new Schwartzchild body with name, m, x, v."},
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
  lua_pushfstring(L, "body (%s)", b->name);
  return 1;
}

static lua_record_metamethod body_mm[] = {
  {"__tostring", body_tostring},
  {NULL, NULL}
};

// This function constructs our "n squared" model. It takes a table of 
// arguments.
static int n2_new(lua_State* L)
{
  if ((lua_gettop(L) != 1) || !lua_istable(L, 1))
    luaL_error(L, "Argument must be a table of parameters and bodies.");

  // Where all our bodies at?
  lua_getfield(L, 1, "bodies");
  if (!lua_istable(L, -1))
    luaL_error(L, "bodies must be a list of body objects.");
  body_array_t* bodies = body_array_new();
  {
    int i = 1;
    bool have_swartzchild_body = false;
    while (true)
    {
      lua_rawgeti(L, -1, i);
      if (lua_isnil(L, -1)) break;
      if (!lua_is_body(L, -1))
        luaL_error(L, "Item %d in bodies is not a body.", i);
      body_t* b = lua_to_body(L, -1);
      if (b->schwartzchild)
      {
        if (have_swartzchild_body)
          luaL_error(L, "Only one body in bodies may be a schwartzchild body.");
        else
          have_swartzchild_body = true;
      }
      body_array_append(bodies, b);
      lua_pop(L, 1);
      ++i;
    }
  }
  lua_pop(L, 1);

  // Are we given a gravitational constant G?
  real_t G = GRAVITATIONAL_CONSTANT;
  lua_getfield(L, 1, "G");
  if (lua_isnumber(L, -1))
  {
    G = lua_to_real(L, -1);
    if (G <= 0.0)
      luaL_error(L, "G must be positive.");
  }

  // Create the thingy and push it to the stack.
  model_t* n2 = n_squared_new(G, bodies);

  // Are we given a fixed dt?
  real_t fixed_dt = 0.0;
  lua_getfield(L, 1, "fixed_dt");
  if (lua_isnumber(L, -1))
  {
    fixed_dt = lua_to_real(L, -1);
    if (fixed_dt <= 0.0)
      luaL_error(L, "fixed_dt must be positive.");
    model_set_max_dt(n2, fixed_dt);
  }

  lua_push_model(L, n2);
  return 1;
}

static int n2_x_probe(lua_State* L)
{
  if (!lua_isstring(L, 1))
    luaL_error(L, "Argument must be the name of a body.");
  const char* name = lua_tostring(L, 1);
  lua_push_probe(L, n_squared_x_probe_new(name));
  return 1;
}

static int n2_v_probe(lua_State* L)
{
  if (!lua_isstring(L, 1))
    luaL_error(L, "Argument must be the name of a body.");
  const char* name = lua_tostring(L, 1);
  lua_push_probe(L, n_squared_v_probe_new(name));
  return 1;
}

static int n2_E_probe(lua_State* L)
{
  lua_push_probe(L, n_squared_E_probe_new());
  return 1;
}

static void lua_register_n_squared(lua_State* L)
{
  // Create a new table and fill it with our gravity models.
  lua_newtable(L);

  lua_pushcfunction(L, n2_new);
  lua_setfield(L, -2, "new");

  lua_pushcfunction(L, n2_x_probe);
  lua_setfield(L, -2, "x_probe");

  lua_pushcfunction(L, n2_v_probe);
  lua_setfield(L, -2, "v_probe");

  lua_pushcfunction(L, n2_E_probe);
  lua_setfield(L, -2, "E_probe");

  lua_setglobal(L, "n_squared");
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

  lua_register_n_squared(L);
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

