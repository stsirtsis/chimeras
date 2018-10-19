CC = gcc
CFLAGS = -O2 -Wall -Wextra
#FCFLAGS += -I/usr/include
PROGRAMS = BLif2DarisSigma_0.2 LIF_2D_Classic LIF_2D_Classic_Parallel2 LIF_2D_Classic_Parallel4 LIF_2D_Classic_Parallel5 Generator

all: $(PROGRAMS)

old: BLif2DarisSigma_0.2

new: LIF_2D_Classic

par2: LIF_2D_Classic_Parallel2
par4: LIF_2D_Classic_Parallel4
par5: LIF_2D_Classic_Parallel5

gen: Generator

LIF_2D_Classic: LIF_2D_Classic.c
		$(CC) $(CFLAGS) LIF_2D_Classic.c -o LIF_2D_Classic

BLif2DarisSigma_0.2: BLif2DarisSigma_0.2.c
		$(CC) $(CFLAGS) BLif2DarisSigma_0.2.c -o BLif2DarisSigma_0.2

LIF_2D_Classic_Parallel2: LIF_2D_Classic_Parallel.c
		$(CC) $(CFLAGS) LIF_2D_Classic_Parallel.c -D NUMTHREADS=2 -pthread -o LIF_2D_Classic_Parallel2

LIF_2D_Classic_Parallel4: LIF_2D_Classic_Parallel.c
		$(CC) $(CFLAGS) LIF_2D_Classic_Parallel.c -D NUMTHREADS=4 -pthread -o LIF_2D_Classic_Parallel4

LIF_2D_Classic_Parallel5: LIF_2D_Classic_Parallel.c
		$(CC) $(CFLAGS) LIF_2D_Classic_Parallel.c -D NUMTHREADS=5 -pthread -o LIF_2D_Classic_Parallel5

Generator: Generator.c
		$(CC) $(CFLAGS) Generator.c -o Generator

clean:
		rm -f *.o *.mod *.MOD
		rm -f *~ $(PROGRAMS) *dat*
		rm -r Results*
		rm initial*
