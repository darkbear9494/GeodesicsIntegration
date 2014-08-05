gcc -o geocross GeoCross.c PseN.c dopri8.c geoIntegrator.c iniParlist.c geoSave.c geoMove.c iniTidal.c geoCollide.c
#gcc -fopenmp -o geocross GeoCross.c PseN.c dopri8.c geoIntegrator.c iniParlist.c geoSave.c geoMove.c iniTidal.c
echo compile succeed!
./geocross
echo remove geocross
rm geocross
