/**
 * @brief matrix calculation
 */
#include <stdlib.h>
#include "../inc/matrix.h"
Matrix mTranspose(Matrix s);

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
            trace = cAdd(trace, ind.pdata[i][i]);
        }
        free(ind.pdata);
        return trace;
    }
    errHandler("trace: incompatible dimensions.");
    return (dComplex){0, 0};
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
void rowScale(Matrix M, int n, dComplex t)
{
    pIndex ind = mapIndex(M);
    for (int i = 0; i < M.col; i++)
    {
        ind.pdata[n][i] = cProd(ind.pdata[n][i], t);
    }
    free(ind.pdata);
}
void rowAddon(Matrix M, int m, int n, dComplex t)
{
    // add t*row n onto row m
    pIndex ind = mapIndex(M);
    for (int i = 0; i < M.col; i++)
    {
        ind.pdata[m][i] = cAdd(ind.pdata[m][i], cProd(ind.pdata[n][i], t));
    }
    free(ind.pdata);
}
void singleRowElim(Matrix M,
                   int m /* row to be eliminated*/,
                   int n /*row refered*/,
                   int p /*col refered*/)
{
    dComplex t;
    pIndex ind = mapIndex(M);
    t.real = -cDiv(ind.pdata[m][p], ind.pdata[n][p]).real;
    t.imag = -cDiv(ind.pdata[m][p], ind.pdata[n][p]).imag;
    rowAddon(M, m, n, t);
    free(ind.pdata);
}
Matrix gaussElim(Matrix M)
{
    dComplex *pCopy = (dComplex *)malloc(M.row * M.col * dSize);
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
                        // debug
                        // printm(M, "%2.1f ", true);
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
                    singleRowElim(M, k, i, j);
                    // debug
                    // printm(M, "%2.1f ", true);
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
            singleRowElim(M, i, pivot[sp].row, pivot[sp].col);
        }
        dComplex scaler = cDiv((dComplex){1, 0}, ind.pdata[pivot[sp].row][pivot[sp].col]);
        rowScale(M, pivot[sp].row, scaler);
        // for (int i = pivot[sp].row; i >= 0; i--)
        // {

        // }
    }
    free(ind.pdata);
    return M;
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
Matrix randMatrix(int row, int col)
{
    dComplex *p = (dComplex *)malloc(dSize * row * col);
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