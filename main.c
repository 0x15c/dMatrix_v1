#include <stdio.h>
#include <math.h>
#include "src/matrix.c"
int main(int argc, char **argv)
{
    dComplex exp = {0, 2 * PI / 3};
    dComplex base = {E, 0};
    dComplex result = cPow(base, exp);
    // printc("%6.4f", result, true);
    dComplex dataBlock[] = {
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {2, 0},
        {0, 0},
        {1, 0},
        {0, 0},
        {3, 0},
    };
    Matrix matrix1 = mInit(3, 3);
    matrix1.pdata = dataBlock;
    printm(matrix1, "%2.1f ", true);
    printf("\n");
    // singleRowElim(matrix1, 1, 0, 0);
    rowElimination(matrix1);
    printm(matrix1, "%2.1f ", true);
    // pIndex matrix1_index = mapIndex(matrix1);
    // dComplex trace = mTrace(matrix1);
    // rowExchange(matrix1, 1, 0);
    // printm(matrix1, "%1.1f ");
    // printf("\n");
    // Matrix product = mProd(matrix1, matrix1);
    // printm(product, "%1.1f ");
    // printf("\n");
    // mConj(product);
    // printm(product, "%1.1f ");
    // printf("\n");
    // printm(mHermitian(product), "%1.1f ");
    return 0;
}