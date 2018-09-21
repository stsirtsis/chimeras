CC = gcc
CFLAGS = -O2 -Wall -Wextra
#FCFLAGS += -I/usr/include
PROGRAMS = BLif2DarisSigma_0.2 LIF_2D_Classic

all: $(PROGRAMS)

old: BLif2DarisSigma_0.2

new: LIF_2D_Classic

%: %.o
		$(CC) $(CFLAGS) -o $@ $^

%.o: %.c
		$(CC) $(CFLAGS) -c $<

clean:
		rm -f *.o *.mod *.MOD
		rm -f *~ $(PROGRAMS) *dat*
