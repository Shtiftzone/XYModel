################################################################################
#
# Makefile to compile and link C programs
#
# Version valid for Linux machines
#
# "make" compiles and links the specified main programs and modules
# using the specified libraries (if any), and produces the executables
#
# "make clean" removes all files created by "make"
#
# Dependencies on included files are automatically taken care of
#
################################################################################

all: rmxeq mkdep mkxeq
.PHONY: all


# Main programs and required modules

MAIN = testlx timelx

MODULES = ranlxs ranlxd ranlux_common


# Additional include directories

INCPATH =


# Additional libraries to be included

LIBS = m

LIBPATH =


# Scheduling and optimization options

CFLAGS = -std=c99 -pedantic -fstrict-aliasing \
         -Wall -Wstrict-prototypes -Werror -O2

# The available program-specific optimization options are -DSSE2 and -DAVX2

############################## do not change ###################################

SHELL=/bin/bash

CC=gcc

PGMS= $(MAIN) $(MODULES)

-include $(addsuffix .d,$(PGMS))


# Rule to make dependencies

$(addsuffix .d,$(PGMS)): %.d: %.c Makefile
	@ $(CC) -MM -ansi $(addprefix -I,$(INCPATH)) $< -o $@


# Rule to compile source programs

$(addsuffix .o,$(PGMS)): %.o: %.c Makefile
	$(CC) $< -c $(CFLAGS) \
        $(addprefix -I,$(INCPATH)) -o $(notdir $@)


# Rule to link object files

$(MAIN): %: %.o $(addsuffix .o,$(MODULES)) Makefile
	$(CC) $< $(addsuffix .o,$(MODULES)) $(CFLAGS) \
        $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBS)) -o $@


# Produce executables

mkxeq: $(MAIN)


# Remove old executables and old error log file

rmxeq:
	@ -rm -f $(MAIN); \
        echo "delete old executables"


# Make dependencies

mkdep:  $(addsuffix .d,$(PGMS))
	@ echo "generate tables of dependencies"


# Clean directory

clean:
	@ -rm -rf *.d *.o $(MAIN)
.PHONY: clean

################################################################################
