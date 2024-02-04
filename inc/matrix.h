#include "complex.h"
void errHandler(const char *errLog);
#define MAX(a, b) ((a) > (b)) ? (a) : (b)
const int dSize = sizeof(dComplex);
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
void printm(Matrix target, const char *format, flag RealOnly);
Matrix mTranspose(Matrix s);
typedef struct Matrix rvec;
Matrix mInit(int row, int col);
void cSwap(dComplex *s, dComplex *t);
void cCpy(dComplex source, dComplex *target);
pIndex mapIndex(Matrix M);
dComplex mTrace(Matrix M);
void rowExchange(Matrix M, int n, int m);
void rowScale(Matrix M, int n, dComplex t);
void rowAddon(Matrix M, int m, int n, dComplex t);
void singleRowElim(Matrix M, int m /* row to be eliminated*/, int n /*row refered*/, int p /*col refered*/);
Matrix gaussElim(Matrix M);
Matrix reducedRowElim(Matrix N);
void colExchange(Matrix M, int n, int m);
Matrix mAdd(Matrix s, Matrix t);
Matrix mProd(Matrix s, Matrix t);
void mConj(Matrix s);
Matrix mTranspose(Matrix s);
Matrix mHermitian(Matrix s);
void printm(Matrix target, const char *format, flag RealOnly);
Matrix randMatrix(int row, int col);
dComplex det(Matrix M);
void errHandler(const char *errLog);