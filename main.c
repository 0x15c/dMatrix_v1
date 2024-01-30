#include <stdio.h>
#include <math.h>
#include "src/complex.c"
int main(int argc, char **argv)
{
    printf("this is a test\n");
    dComplex a = {2.0, 2.0};
    dComplex k = cAdd(a, (dComplex){1.0, 2.0});
    dComplex d = cConj(a);
    printc("%3.2f", k, true);
    dComplex mw = cDiv((dComplex){1.0, 2.0}, (dComplex){2.0, 1.0});
    dEuler e = complex2Euler(a);
    return 0;
}