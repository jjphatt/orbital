# Makefile -- Use this to build on *NIX systems.

# Options set on command line.
verbose    = not-set
polymec    = not-set
ide        = not-set

# This proxies everything to the builddir cmake.

cputype = $(shell uname -m | sed "s/\\ /_/g")
systype = $(shell uname -s)

BUILDDIR := build/$(systype)-$(cputype)
CONFIG_FLAGS = -DUNIX=1

# Process configuration options.

# Verbose builds?
ifeq ($(verbose), 1)
  CONFIG_FLAGS += -DCMAKE_VERBOSE_MAKEFILE=1
endif

BUILDDIR := ${BUILDDIR}-`basename ${CC}`

# Location of polymec library.
ifneq ($(polymec), not-set)
  CONFIG_FLAGS += -DPOLYMEC_PREFIX:PATH=$(polymec)
else
  CONFIG_FLAGS += -DPOLYMEC_PREFIX:PATH=/usr/local
endif

# Integrated Development Environment (IDE)?
ifeq ($(ide), codeblocks)
  CONFIG_FLAGS += -G "CodeBlocks - Unix Makefiles"
else
  ifeq ($(ide), xcode)
    CONFIG_FLAGS += -GXcode
  else
    ifeq ($(ide), codelite)
      CONFIG_FLAGS += -G "CodeLite - Unix Makefiles"
    else
      ifeq ($(ide), eclipse)
        CONFIG_FLAGS += -G "Eclipse CDT4 - Unix Makefiles"
      else
        ifeq ($(ide), kate)
          CONFIG_FLAGS += -GKate 
        else
          ifeq ($(ide), sublime)
            CONFIG_FLAGS += -G "Sublime Text 2 - Unix Makefiles" 
          endif
        endif
      endif
    endif
  endif
endif

define run-config
@mkdir -p $(BUILDDIR)
@cd $(BUILDDIR) && cmake $(CURDIR) $(CONFIG_FLAGS)
endef

all test clean install:
	@if [ ! -f $(BUILDDIR)/Makefile ]; then \
		more INSTALL; \
	else \
		make -C $(BUILDDIR) $@ $(MAKEFLAGS); \
	fi

doc:
	make -C doc $@ $(MAKEFLAGS); \

config: distclean
	$(run-config)

distclean:
	@rm -rf $(BUILDDIR)

.PHONY: config distclean all clean install uninstall 
