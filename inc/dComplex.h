#ifndef DCOMPLEX_H
#define DCOMPLEX_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define PI 3.141592653589793238462643383279502884L
#define E 2.718281828459045235360287471352662498L
#define EPS 1e-6
typedef enum posQuadrant
{
    Quadrant_1_4,
    Quadrant_2,
    Quadrant_3,
    RealAxis,
    ImagAxis,
    Origin,
} cPos;
typedef enum flagTruth
{
    false = 0,
    true = 1,
} flag;
typedef struct cartesianComplex
{
    float real, imag;
} dComplex;
typedef struct eulerComplex
{
    float modulus;
    float argument;
} dEuler;
cPos getPos(dComplex s);
void printc(const char *format, dComplex c, flag RealOnly);
dComplex cAdd(dComplex s1, dComplex s2);
dComplex cMin(dComplex s1, dComplex s2);
dComplex cConj(dComplex s);
dComplex cProd(dComplex s1, dComplex s2);
float cModu(dComplex s);
dComplex cScale(dComplex s, float factor);
dComplex cDiv(dComplex s1, dComplex s2);
dEuler complex2Euler(dComplex s);
dComplex euler2Complex(dEuler s);
dComplex cPow(dComplex base, dComplex exp);
#endif