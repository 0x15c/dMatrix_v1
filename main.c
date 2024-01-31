#include <stdio.h>
#include <math.h>
#include "src/matrix.c"
int main(int argc, char **argv)
{
    dComplex exp = {0, 2 * PI / 3};
    dComplex base = {E, 0};
    dComplex result = cPow(base, exp);
    printc("%6.4f", result, true);
    dComplex dataBlock[] = {
        {1, 1},
        {1, 2},
        {1, 3},
        {2, 1},
        {2, 2},
        {2, 3},
        {3, 1},
        {3, 2},
        {3, 3},
    };
    Matrix matrix1 = mInit(3, 3);
    matrix1.pdata = dataBlock;
    printm(matrix1, "%1.1f ");
    printf("\n");

    pIndex matrix1_index = mapIndex(matrix1);
    dComplex trace = mTrace(matrix1);
    rowExchange(matrix1, 1, 0);
    printm(matrix1, "%1.1f ");
    printf("\n");
    Matrix product = mProd(matrix1, matrix1);
    printm(product, "%1.1f ");
    printf("\n");
    mConj(product);
    printm(product, "%1.1f ");
    printf("\n");
    printm(mHermitian(product), "%1.1f ");
    return 0;
}