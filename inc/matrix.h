#ifndef DCOMPLEX_H
#include "dComplex.h"
#endif
#ifndef MATRIX_H
#define MATRIX_H
void errHandler(const char *errLog);
#define MAX(a, b) ((a) > (b)) ? (a) : (b)
#define dComplexSize sizeof(dComplex)
typedef struct Matrix
{
    const int row, col;
    dComplex *pdata;
} Matrix;
typedef struct entryIndex
{
    int row;
    int col;
} entryIndex;
typedef struct ptrIndex
{
    const int row, col;
    dComplex **pdata;
} pIndex;
int mRank(Matrix N);
flag mOperationMode(Matrix M /** @return 0, real only; 1, complex*/);
Matrix mTranspose(Matrix s);
typedef struct Matrix rvec;
Matrix mInit(int row, int col);
void cSwap(dComplex *s, dComplex *t, const flag isComplex);
void cCpy(dComplex source, dComplex *target);
pIndex mapIndex(Matrix M);
dComplex mTrace(Matrix M);
void rowExchange(Matrix M, int n, int m, const flag isComplex);
void rowScale(Matrix M, int n, dComplex t, const flag isComplex);
void rowAddon(Matrix M, int m, int n, dComplex t, const flag isComplex);
void singleRowElim(Matrix M, int m /* row to be eliminated*/, int n /*row refered*/, int p /*col refered*/, const flag isComplex);
Matrix gaussElim(Matrix M);
Matrix reducedRowElim(Matrix N);
// void colExchange(Matrix M, int n, int m);
Matrix mAdd(Matrix s, Matrix t);
Matrix mProd(Matrix s, Matrix t);
void mConj(Matrix s);
Matrix mTranspose(Matrix s);
Matrix mHermitian(Matrix s);
void printm(Matrix target, const char *format);
Matrix randMatrix(int row, int col);
dComplex det(Matrix M);
void errHandler(const char *errLog);
Matrix mZeros(int row, int col);
Matrix mSetVal(int row, int col, dComplex val);
Matrix mIdentity(int dim);
void printm_pretty(Matrix M);
#endif