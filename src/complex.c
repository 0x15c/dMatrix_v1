/**
 * @brief complex arithmetic support
 *
 *
 *
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI 3.141592653589793238462643383279502884L
#define E 2.718281828459045235360287471352662498L
typedef enum
{
    Quadrant_1_4,
    Quadrant_2,
    Quadrant_3,
    RealAxis,
    ImagAxis,
    Origin,
} cPos;
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
cPos getPos(dComplex s)
{
    if (s.real > 0 && s.imag != 0)
        return Quadrant_1_4;
    else if (s.real < 0 && s.imag > 0)
        return Quadrant_2;
    else if (s.real < 0 && s.imag < 0)
        return Quadrant_3;
    else if (s.real == 0 && s.imag != 0)
        return ImagAxis;
    else if (s.imag == 0 && s.real != 0)
        return RealAxis;
    else
        return Origin;
}
void printc(const char *format, dComplex c, flag isCR)
{
    printf(format, c.real);
    printf(" + j");
    printf(format, c.imag);
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
    switch (getPos(s))
    {
    case Quadrant_1_4:
        // Quadrant I and IV
        return (dEuler){
            cModu(s),
            atan(s.imag / s.real),
        };
    case Quadrant_2:
        // Quadrant II
        return (dEuler){
            cModu(s),
            PI + atan(s.imag / s.real),
        };
    case Quadrant_3:
        // Quadrant III
        return (dEuler){
            cModu(s),
            -PI + atan(s.imag / s.real),
        };
    case ImagAxis:
        return (dEuler){
            fabs(s.imag),
            s.imag > 0 ? PI / 2 : -PI / 2,
        };
    case RealAxis:
        return (dEuler){
            fabs(s.real),
            s.real > 0 ? 0 : PI,
        };
    case Origin:
        return (dEuler){0, 0};
    }
}
dComplex euler2Complex(dEuler s)
{
    return (dComplex){
        s.modulus * cos(s.argument),
        s.modulus * sin(s.argument),
    };
}
dComplex cPow(dComplex base, dComplex exp)
{
    /**
     * @return gives the value of base^exp, in principle argument interval
     * @note this function may give unwanted value because of the multi-value nature
     */
    dEuler temp;
    temp.modulus = powf(complex2Euler(base).modulus, exp.real) * powf(E, -(complex2Euler(base).argument));
    temp.argument = exp.imag * log(complex2Euler(base).modulus) + exp.real * complex2Euler(base).argument;
    return euler2Complex(temp);
}