FC = gfortran
FCFLAGS = -O2
FCFLAGS += -I/usr/include
PROGRAMS = lif2d-classic-new lif2d-classic

all: $(PROGRAMS)

old: lif2d-classic

new: lif2d-classic-new

%: %.o
		$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
		$(FC) $(FCFLAGS) -c $<

%.o: %.f
		$(FC) $(FCFLAGS) -c $<

clean:
		rm -f *.o *.mod *.MOD
		rm -f *~ $(PROGRAMS) *dat*
