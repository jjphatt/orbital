// Copyright (c) 2016, Jeffrey N. Johnson
// All rights reserved.

#include "orbital.h"
#include "lua_orbital.h"

int main(int argc, char* argv[]) 
{
  return lua_driver(argc, argv, lua_register_orbital_modules);
}

