/**
 * @brief complex arithmetic support
 *
 *
 *
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
typedef enum
{
    false = 0,
    true = 1,
} flag;
typedef struct
{
    float real, imag;
} dComplex;
typedef struct
{
    float modulus;
    float argument;
} dEuler;
void printc(const char *printArg, dComplex c, flag isCR)
{
    printf(printArg, c.real);
    printf(" + j");
    printf(printArg, c.imag);
    if (isCR)
        printf("\n");
}

dComplex cAdd(dComplex s1, dComplex s2)
{
    s1.real += s2.real;
    s1.imag += s2.imag;
    return s1;
}

dComplex cMin(dComplex s1, dComplex s2)
{
    s1.real -= s2.real;
    s1.imag -= s2.imag;
    return s1;
}

dComplex cConj(dComplex s)
{
    s.imag = -s.imag;
    return s;
}

dComplex cProd(dComplex s1, dComplex s2)
{
    return (dComplex){
        s1.real * s2.real - s1.imag * s2.imag,
        s1.real * s2.imag + s1.imag * s2.real,
    };
}

float cModu(dComplex s)
{
    return sqrt(cProd(s, cConj(s)).real);
}

dComplex cScale(dComplex s, float factor)
{
    return (dComplex){
        s.real * factor,
        s.imag * factor,
    };
}

dComplex cDiv(dComplex s1, dComplex s2)
{
    return (cProd(s1, (cScale(cConj(s2), 1 / cProd(s2, cConj(s2)).real))));
}
dEuler complex2Euler(dComplex s)
{
    return (dEuler){
        cModu(s),
        atan(s.imag / s.real),
    };
}