objects = GeoCross.o dopri8.o PseN.o geoIntegrator.o iniPar.o
edit: $(objects)
	cc -o edit $(objects)

GeoCross.o: GeoCross.h dopri8.h PseN.h
geoIntegrator.o: GeoCross.h dopri8.h PseN.h
iniPar.o: GeoCross.h
PseN.o: PseN.h
dopri8.o: dopri8.h

.PHONY: clean
clean:
	rm edit $(objects)
