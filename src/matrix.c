/**
 * @brief matrix calculation
 */
#include <stdlib.h>
#include "complex.c"
void errHandler(const char *errLog);
const int dSize = sizeof(dComplex);
typedef struct
{
    const int row, col;
    dComplex *pdata;
} Matrix;
typedef struct
{
    const int row, col;
    dComplex **pdata;
} pIndex;

typedef struct Matrix rvec;
Matrix mInit(int row, int col)
{
    // set aside memory for new created matrix
    dComplex *pdata = (dComplex *)malloc(row * col * dSize);
    return (Matrix){
        .row = row,
        .col = col,
        .pdata = pdata,
    };
}
pIndex mapIndex(Matrix M)
{
    dComplex **p = (dComplex **)malloc(M.row * sizeof(dComplex *));
    for (int i = 0; i < M.row; ++i)
    {
        p[i] = M.pdata + i * M.col;
    }
    return (pIndex){
        M.row,
        M.col,
        p,
    };
}
dComplex mTrace(Matrix M)
{
    if (M.col == M.row)
    {
        pIndex ind = mapIndex(M);
        dComplex trace = {0, 0};
        for (int i = 0; i < M.row; i++)
        {
            trace.real += ind.pdata[i][i].real;
            trace.imag += ind.pdata[i][i].imag;
        }
        free(ind.pdata);
        return trace;
    }
    errHandler("trace: incompatible dimensions, abort.");
}

void errHandler(const char *errLog)
{
    printf("%s\n", errLog);
}