# Makefile for generic compilation -- includes a hack
# to get the dependencies from .h files directly from the compiler;
#     GNU compilers: AUTO -> -MM
# non-GNU compilers: AUTO -> -M
# ldd/00

MAKEFILE = Makefile

CC = g++

# additional include directories
INCPATH =

# additional libraries to be included  
LIBS = 
LIBPATH =

# scheduling and optimization options  
#CFLAGS = -O3 -march=pentium4 -mfpmath=sse -msse2 -malign-double
CFLAGS = -O2 

#macros (such as -DSSE -DSSE2 -DP4)
MACRO = -DSSE2 -DAVX2

AUTO = -MM

SRCS = main.cpp ranlxs.cpp MultiHistRW.cpp Autocorr.cpp

OBJS = $(SRCS:.cpp=.o)

default: Rw 

%.d: %.cpp
	set -e; $(CC) $(AUTO) $< \
                  | sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
                [ -s $@ ] || rm -f $@

%.o: %.cpp $(MAKEFILE) 
	$(CC) -c $(CFLAGS) $(MACRO) $< -o $@

include $(SRCS:.cpp=.d)

.PHONY: clean

Rw: $(MAKEFILE) $(OBJS)
	$(CC) $(CFLAGS) $(MACRO) -o $@ $(OBJS)

Autoc: $(MAKEFILE) $(OBJS)
	$(CC) $(CFLAGS) $(MACRO) -o $@ $(OBJS)

clean:
	rm -f core *.o *.d *~

tags:
	find . -name "*.cpp" -or -name "*.h" | etags -
