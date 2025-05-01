#
#  Makefile
#
#  XY model (2D) with Wolff updates and filtration-based topological analysis
#

CFLAGS   = -std=c11 -O3 -Wall -pedantic -Iranlux-3.4
RANLUX   = ranlux-3.4/ranlxs.o ranlux-3.4/ranlxd.o ranlux-3.4/ranlux_common.o
LIBS     = -lm

L_VALUES = 64 96 128 192 256 384 512 768 1024

all: $(L_VALUES)

$(L_VALUES):
	@echo "Compiling for L=$@"
	gcc $(CFLAGS) -DDIM=2 -DCELL -DL=$@ xy.c $(RANLUX) $(LIBS) -o xy2d-$@

clean:
	rm -f *.o xy2d-*
