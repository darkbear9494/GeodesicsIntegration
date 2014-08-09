gcc -o geocross GeoCross.c PseN.c dopri8.c geoIntegrator.c iniParlist.c geoSave.c geoMove.c iniTidal.c geoCollide.c geoToy.c iniQuant.c
#gcc -fopenmp -o geocross GeoCross.c PseN.c dopri8.c geoIntegrator.c iniParlist.c geoSave.c geoMove.c iniTidal.c
echo geocross: COMPILE SUCCEED!
./geocross
#echo remove geocross
rm geocross
