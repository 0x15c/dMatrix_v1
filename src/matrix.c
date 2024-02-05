/**
 * @brief matrix calculation
 */
#include <stdlib.h>
#include "matrix.h"
flag mOperationMode(Matrix M /** @return 0, real only; 1, complex*/)
{
    float acc = 0.0f;
    for (int i = 0; i < M.col * M.row; i++)
    {
        acc += fabs(M.pdata[i].imag);
    }
    if (acc < EPS)
    {
        return 0;
    }
    return 1;
}
Matrix mInit(int row, int col)
{
    // set aside memory for new created matrix
    dComplex *pdata = (dComplex *)malloc(row * col * dComplexSize);
    return (Matrix){row, col, pdata};
}
void cSwap(dComplex *s, dComplex *t, const flag isComplex)
{
    dComplex temp;
    temp.real = s->real;
    s->real = t->real;
    t->real = temp.real;
    if (isComplex)
    {
        temp.imag = s->imag;
        s->imag = t->imag;
        t->imag = temp.imag;
    }
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
            trace = cAdd(trace, ind.pdata[i][i]);
        }
        free(ind.pdata);
        return trace;
    }
    errHandler("trace: incompatible dimensions.");
    return (dComplex){0, 0};
}
void rowExchange(Matrix M, int n, int m, const flag isComplex)
{
    pIndex ind = mapIndex(M);
    for (int i = 0; i < M.col; i++)
    {
        cSwap(&ind.pdata[n][i], &ind.pdata[m][i], isComplex);
    }
    free(ind.pdata);
}
void rowScale(Matrix M, int n, dComplex t, const flag isComplex)
{
    pIndex ind = mapIndex(M);
    for (int i = 0; i < M.col; i++)
    {
        ind.pdata[n][i] = isComplex ? cProd(ind.pdata[n][i], t)
                                    : (dComplex){t.real * ind.pdata[n][i].real, 0};
    }
    free(ind.pdata);
}
void rowAddon(Matrix M, int m, int n, dComplex t, const flag isComplex)
{
    // add t*row n onto row m
    pIndex ind = mapIndex(M);
    for (int i = 0; i < M.col; i++)
    {
        ind.pdata[m][i] = isComplex ? cAdd(ind.pdata[m][i],
                                           cProd(ind.pdata[n][i], t))
                                    : cAdd(ind.pdata[m][i],
                                           (dComplex){ind.pdata[n][i].real * t.real, 0});
    }
    free(ind.pdata);
}
void singleRowElim(Matrix M,
                   int m /* row to be eliminated*/,
                   int n /*row refered*/,
                   int p /*col refered*/,
                   const flag isComplex)
{
    dComplex t;
    pIndex ind = mapIndex(M);
    t.real = -cDiv(ind.pdata[m][p], ind.pdata[n][p]).real;
    t.imag = -cDiv(ind.pdata[m][p], ind.pdata[n][p]).imag;
    rowAddon(M, m, n, t, isComplex);
    free(ind.pdata);
}
Matrix gaussElim(Matrix M)
{
    const flag isComplex = mOperationMode(M);
    dComplex *pCopy = (dComplex *)malloc(M.row * M.col * dComplexSize);
    for (int i = 0; i < M.row * M.col; i++)
        cCpy(M.pdata[i], &pCopy[i]);
    ;
    dComplex *pTemp = M.pdata;
    M.pdata = pCopy;
    pIndex ind = mapIndex(M);
    int i, j, k;
    i = 0;
    for (j = i; j < M.col && i < M.row;)
    {
        // starting with current pivot position (if any)
        for (; i < M.row; i++)
        {
            // row index increment
            if ((isComplex ? cModu(ind.pdata[i][j]) : fabs(ind.pdata[i][j].real)) < EPS)
            {
                /**
                 * if encounters zero element, search for next non-zero element and exchange recursively
                 */
                for (k = i + 1; k < M.row; k++)
                {
                    if ((isComplex ? cModu(ind.pdata[k][j]) : fabs(ind.pdata[k][j].real)) > EPS)
                    {
                        rowExchange(M, i, k, isComplex);
                        // debug
                        // printm(M, "%2.1f ");
                        // printf("\n");
                        break;
                    }
                    // break;
                }
            }
            if (k == M.row - 1)
            {
                j++;
                break;
                // encounters zero column
            }
            else
            {
                /**
                 * now it is gurenteed to have rest elements with non-zero entry in current column
                 * safe to do the single row elimination
                 */
                for (k = i + 1; k < M.row; k++)
                {
                    singleRowElim(M, k, i, j, isComplex);
                    // debug
                    // printm(M, "%2.1f ");
                    // printf("\n");
                }
                j++;
            }
        }
    }
    M.pdata = pTemp;
    return (Matrix){M.row, M.col, pCopy};
}
Matrix reducedRowElim(Matrix N)
{
    const flag isComplex = mOperationMode(N);
    Matrix M = gaussElim(N);
    entryIndex pivot[MAX(M.row, M.col)];
    pIndex ind = mapIndex(M);
    int sp = 0;
    int rowSearched, colSearched;
    rowSearched = colSearched = 0;
    for (int j = colSearched; j < M.col; j++)
    {
        for (int i = rowSearched; i < M.row; i++)
        {
            if (cModu(ind.pdata[i][j]) > EPS) // hit a pivot
            {
                pivot[sp].row = i;
                pivot[sp].col = j;
                sp++;
                rowSearched++;
                colSearched++;
                break;
            }
        }
    }
    // obtained pivot position indicator, then perform back elimination
    while (sp > 0)
    {
        sp--;
        for (int i = pivot[sp].row - 1; i >= 0; i--)
        {
            singleRowElim(M, i, pivot[sp].row, pivot[sp].col, isComplex);
        }
        dComplex scaler = isComplex ? cDiv((dComplex){1, 0}, ind.pdata[pivot[sp].row][pivot[sp].col])
                                    : (dComplex){1 / ind.pdata[pivot[sp].row][pivot[sp].col].real, 0};
        rowScale(M, pivot[sp].row, scaler, isComplex);
        // for (int i = pivot[sp].row; i >= 0; i--)
        // {

        // }
    }
    free(ind.pdata);
    return M;
}
// void colExchange(Matrix M, int n, int m)
// {
//     pIndex ind = mapIndex(M);
//     for (int i = 0; i < M.col; i++)
//     {
//         cSwap(&ind.pdata[i][n], &ind.pdata[i][m]);
//     }
//     free(ind.pdata);
// }
Matrix mAdd(Matrix s, Matrix t)
{
    if (s.row == t.row && s.col == t.col)
    {
        dComplex *r = (dComplex *)malloc(dComplexSize * s.row * t.col);

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
    flag isComplex = mOperationMode(s);
    if (s.col == t.row)
    {
        dComplex *r = (dComplex *)malloc(dComplexSize * s.row * t.col);
        for (int i = 0; i < s.row; i++)
        {
            for (int j = 0; j < t.col; j++)
            {
                r[i * s.row + j].imag = r[i * s.row + j].real = 0;
                for (int k = 0; k < s.col; k++)
                {
                    r[i * s.row + j] = cAdd(r[i * s.row + j],
                                            isComplex ? cProd(s.pdata[i * s.row + k],
                                                              t.pdata[k * t.col + j])
                                                      : (dComplex){
                                                            s.pdata[i * s.row + k].real * t.pdata[k * t.col + j].real,
                                                            0});
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
    dComplex *p = (dComplex *)malloc(dComplexSize * s.row * s.col);
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
void printm(Matrix target, const char *format)
{
    int i, j = 0;
    flag realOnly = !mOperationMode(target);
    for (i = 0; i < target.row; i++)
    {
        for (j = 0; j < target.col; j++)
        {
            printc(format, target.pdata[j + i * target.col], realOnly);
        }
        printf("\n");
    }
}
Matrix randMatrix(int row, int col)
{
    dComplex *p = (dComplex *)malloc(dComplexSize * row * col);
    for (int i = 0; i < row * col; i++)
    {
        p[i].real = (rand() % 16);
        p[i].imag = (rand() % 16);
    }
    return (Matrix){row, col, p};
}
dComplex det(Matrix M)
{
    Matrix N = gaussElim(M);
    pIndex ind = mapIndex(N);
    dComplex det = {1, 0};
    for (int i = 0; i < N.row; i++)
    {
        det = cProd(det, ind.pdata[i][i]);
    }
    free(ind.pdata);
    free(N.pdata);
    return det;
}
void errHandler(const char *errLog)
{
    printf("%s\n", errLog);
    // exit(1);
}