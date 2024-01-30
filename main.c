#include <stdio.h>
#include <math.h>
#include "src/complex.c"
int main(int argc, char **argv)
{
    dComplex exp = {0, 2 * PI / 3};
    dComplex base = {E, 0};
    dComplex result = cPow(base, exp);
    printc("%6.3f", result, true);
    return 0;
}