FC = gfortran
FCFLAGS = -O2
FCFLAGS += -I/usr/include
PROGRAMS = lif2d-classic-new

all: $(PROGRAMS)

%: %.o
		$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
		$(FC) $(FCFLAGS) -c $<

clean:
		rm -f *.o *.mod *.MOD
		rm -f *~ $(PROGRAMS) *minusp00
