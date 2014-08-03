gcc -o geocross GeoCross.c PseN.c dopri8.c geoIntegrator.c iniPar.c \
geoSave.c geoMove.c
echo compile succeed!
./geocross
echo remove geocross
rm geocross
