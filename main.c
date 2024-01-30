#include <stdio.h>
#include <math.h>
#include "src/complex.c"
int main(int argc, char **argv)
{
    dComplex b = {2, 1};
    dComplex c = cPow(b, b);
    return 0;
}