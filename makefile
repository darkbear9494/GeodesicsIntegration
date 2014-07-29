objects = GeoInte.o dopri8.o PseN.o
edit: $(objects)
	cc -o edit $(objects)

GeoInte.o: Para.h dopri8.h PseN.h
PseN.o: PseN.h Para.h
dopri8.o: dopri8.h

.PHONY: clean
clean:
	rm edit $(objects)
