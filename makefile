objects = GeoInte.o odeint.o bsstep.o mmid.o pzextr.o
edit: $(objects)
	cc -o edit $(objects)

GeoInte.o: GeoInte.h
odeint.o bsstep.o mmid.o pzextr.o: nrutil.h

.PHONY: clean
clean:
	rm edit $(objects)
