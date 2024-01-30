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
void cSwap(dComplex *s, dComplex *t)
{
    dComplex temp;
    temp.real = s->real;
    temp.imag = s->imag;
    s->real = t->real;
    s->imag = t->imag;
    t->real = temp.real;
    t->imag = temp.imag;
}
void cCpy(dComplex source, dComplex *target)
{
    target->real = source.real;
    target->imag = source.imag;
}
pIndex mapIndex(Matrix M)
{
    dComplex **p = (dComplex **)malloc(M.row * sizeof(dComplex *)); // using dynamic memory allocation
    // dComplex *p[M.row];
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
    errHandler("trace: incompatible dimensions.");
}
void rowExchange(Matrix M, int n, int m)
{
    pIndex ind = mapIndex(M);
    for (int i = 0; i < M.col; i++)
    {
        cSwap(&ind.pdata[n][i], &ind.pdata[m][i]);
    }
    free(ind.pdata);
}
void colExchange(Matrix M, int n, int m)
{
    pIndex ind = mapIndex(M);
    for (int i = 0; i < M.col; i++)
    {
        cSwap(&ind.pdata[i][n], &ind.pdata[i][m]);
    }
    free(ind.pdata);
}
Matrix mAdd(Matrix s, Matrix t)
{
    if (s.row == t.row && s.col == t.col)
    {
        dComplex *r = (dComplex *)malloc(dSize * s.row * t.col);

        for (int i = 0; i < s.col * s.row; i++)
        {
            cCpy(cAdd(s.pdata[i], t.pdata[i]), &r[i]);
        }
        return (Matrix){
            s.row,
            s.col,
            r,
        };
    }
    errHandler("mAdd: operand dimension unmatch.");
}
Matrix mProd(Matrix s, Matrix t)
{
    if (s.col == t.row)
    {
        dComplex *r = (dComplex *)malloc(dSize * s.row * t.col);
        for (int i = 0; i < s.row; i++)
        {
            for (int j = 0; j < t.col; j++)
            {
                r[i * s.row + j].imag = r[i * s.row + j].real = 0;
                for (int k = 0; k < s.col; k++)
                {
                    r[i * s.row + j] = cAdd(r[i * s.row + j],
                                            cProd(s.pdata[i * s.row + k],
                                                  t.pdata[k * t.col + j]));
                }
            }
        }
        return (Matrix){
            s.row,
            t.col,
            r,
        };
    }
    errHandler("mProd: operand dimension unmatch.");
}
void printm(Matrix target, const char *format)
{
    int i, j = 0;
    for (i = 0; i < target.row; i++)
    {
        for (j = 0; j < target.col; j++)
        {
            printc(format, target.pdata[j + i * target.col], false);
            // printf("[%d,%d]=%f ", i, j, tMat->mat_index[j + i * tMat->col]);
        }
        printf("\n");
    }
}
void errHandler(const char *errLog)
{
    printf("%s\n", errLog);
    exit(1);
}