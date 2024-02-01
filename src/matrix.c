/**
 * @brief matrix calculation
 */
#include <stdlib.h>
#include "complex.c"

void errHandler(const char *errLog);
const int dSize = sizeof(dComplex);
typedef struct Matrix
{
    const int row, col;
    dComplex *pdata;
} Matrix;
void printm(Matrix target, const char *format, flag RealOnly);
typedef struct ptrIndex
{
    const int row, col;
    dComplex **pdata;
} pIndex;
typedef struct Matrix rvec;

Matrix mInit(int row, int col)
{
    // set aside memory for new created matrix
    dComplex *pdata = (dComplex *)malloc(row * col * dSize);
    return (Matrix){row, col, pdata};
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
    return (pIndex){M.row, M.col, p};
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
    return (dComplex){0, 0};
}
void rowExchange(Matrix M, int n, int m)
{
    for (int i = 0; i < M.col; i++)
    {
        cSwap(&M.pdata[n * M.row + i], &M.pdata[m * M.row + i]);
    }
}
void rowScale(Matrix M, int n, dComplex t)
{
    for (int i = 0; i < M.col; i++)
    {
        M.pdata[n * M.row + i] = cProd(M.pdata[n * M.row + i], t);
    }
}
void rowAddon(Matrix M, int m, int n, dComplex t)
{
    // add t*row n onto row m
    for (int i = 0; i < M.col; i++)
    {
        M.pdata[m * M.row + i] = cAdd(M.pdata[m * M.row + i], cProd(M.pdata[n * M.row + i], t));
    }
}
void singleRowElim(Matrix M,
                   int m /* row to be eliminated*/,
                   int n /*row refered*/,
                   int p /*col refered*/)
{
    dComplex t;
    t.real = -cDiv(M.pdata[m * M.row + p], M.pdata[n * M.row + p]).real;
    t.imag = -cDiv(M.pdata[m * M.row + p], M.pdata[n * M.row + p]).imag;
    rowAddon(M, m, n, t);
}
void rowElimination(Matrix M)
{
    pIndex ind = mapIndex(M);
    int i, j, k;
    i = 0;
    for (j = i; j < M.col;)
    {
        // starting with current pivot position (if any)
        for (; i < M.row; i++)
        {
            // row index increment
            if (cModu(ind.pdata[i][j]) < EPS)
            {
                /**
                 * if encounters zero element, search for next non-zero element and exchange recursively
                 * this is kind of like bubble sort
                 */
                for (k = i + 1; k < M.row; k++)
                {
                    if (cModu(ind.pdata[k][j]) > EPS)
                    {
                        rowExchange(M, i, k);
                        // printm(M, "%2.1f ", true);
                        // printf("\n");
                    }
                }
            }
            if (k == M.row - 1)
                break;
            j++; // encounters zero column
            break;
            /**
             * now it is gurenteed to have rest elements with non-zero entry in current column
             * safe to do the single row elimination
             */
            for (k = i + 1; k < M.row; k++)
            {
                singleRowElim(M, k, i, j);
                // printm(M, "%2.1f ", true);
                // printf("\n");
            }
            j++;
        }
    }
}
// int pivotSearch(int)
// {
//     ;
// }
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
        return (Matrix){s.row, s.col, r};
    }
    errHandler("mAdd: operand dimension unmatch.");
    return (Matrix){0, 0, NULL};
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
        return (Matrix){s.row, t.col, r};
    }
    errHandler("mProd: operand dimension unmatch.");
    return (Matrix){0, 0, NULL};
}
void mConj(Matrix s)
{
    for (int i = 0; i < s.col * s.row; i++)
    {
        s.pdata[i] = cConj(s.pdata[i]);
    }
}
Matrix mTranspose(Matrix s)
{
    dComplex *p = (dComplex *)malloc(dSize * s.row * s.col);
    for (int i = 0; i < s.row; i++)
    {
        for (int j = 0; j < s.col; j++)
        {
            p[i + j * s.row] = s.pdata[j + i * s.col];
        }
    }
    return (Matrix){s.col, s.row, p};
}
Matrix mHermitian(Matrix s)
{
    mConj(s);
    return mTranspose(s);
}
void printm(Matrix target, const char *format, flag RealOnly)
{
    int i, j = 0;
    for (i = 0; i < target.row; i++)
    {
        for (j = 0; j < target.col; j++)
        {
            printc(format, target.pdata[j + i * target.col], RealOnly);
        }
        printf("\n");
    }
}
void errHandler(const char *errLog)
{
    printf("%s\n", errLog);
    // exit(1);
}