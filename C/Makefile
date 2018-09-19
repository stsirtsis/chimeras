CC = gcc
CFLAGS = -O2
#FCFLAGS += -I/usr/include
PROGRAMS = BLif2DarisSigma_0.2 BLif2DarisSigma_0.2_new

all: $(PROGRAMS)

old: BLif2DarisSigma_0.2

new: BLif2DarisSigma_0.2_new

%: %.o
		$(CC) $(CFLAGS) -o $@ $^ 

%.o: %.c
		$(CC) $(CFLAGS) -c $<

clean:
		rm -f *.o *.mod *.MOD
		rm -f *~ $(PROGRAMS) *dat*
