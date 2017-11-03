// Copyright (c) 2016, Jeffrey N. Johnson
// All rights reserved.

#ifndef LUA_ORBITAL_H
#define LUA_ORBITAL_H

#include "lua.h"
#include "orbital.h"

// Registers a complete set of Lua modules for orbital with the given Lua interpreter.
int lua_register_orbital_modules(lua_State* L);

// Pushes a body onto the stack of the given Lua interpreter.
void lua_push_body(lua_State* L, body_t* b);

// Returns true if the variable at the given index in the interpreter L is a 
// body, false otherwise.
bool lua_is_body(lua_State* L, int index);

// Returns a pointer to a body at the given index in L, or NULL if 
// the item there is not a material.
body_t* lua_to_body(lua_State* L, int index);

#endif
